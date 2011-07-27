//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CROSS_STENCIL_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CROSS_STENCIL_HPP

#include <cassert>

#include <algorithm>
#include <numeric>

#include <boost/concept/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/LOOP.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>

namespace PhysBAM
{

template< class T, int D >
struct CROSS_STENCIL
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef T SCALAR_TYPE;

    typedef T VALUE_TYPE;
    static const int DIMENSION = D;
    static const int N_VALUES = 1 + 2*D;

    T values[N_VALUES];

    static CROSS_STENCIL Construct_Zero();

    void Zero();
    int N_Zero() const;
    int N_Nonzero() const;

    T Sum() const;

    //CROSS_STENCIL& operator=(const CROSS_STENCIL& other);
    CROSS_STENCIL& operator+=(const CROSS_STENCIL& other);

    void Add(const VECTOR<T,D>& beta_dv_over_dx_dx, const MULTI_INDEX_TYPE& cell_multi_offset);

    template< class T_STENCIL_PROXY >
    typename boost::enable_if<
        boost::is_same< typename T_STENCIL_PROXY::STENCIL_TYPE, CROSS_STENCIL >
    >::type Assign(const MULTI_INDEX_TYPE& base_multi_index, const T_STENCIL_PROXY& stencil_proxy);
    template< class T_STENCIL_PROXY >
    typename boost::disable_if<
        boost::is_same< typename T_STENCIL_PROXY::STENCIL_TYPE, CROSS_STENCIL >
    >::type Assign(const MULTI_INDEX_TYPE base_multi_index, const T_STENCIL_PROXY stencil_proxy);

    template< class T_STENCIL_PROXY >
    typename boost::enable_if<
        boost::is_same< typename T_STENCIL_PROXY::STENCIL_TYPE, CROSS_STENCIL >
    >::type Add(const MULTI_INDEX_TYPE& base_multi_index, const T_STENCIL_PROXY& stencil_proxy);
    template< class T_STENCIL_PROXY >
    typename boost::disable_if<
        boost::is_same< typename T_STENCIL_PROXY::STENCIL_TYPE, CROSS_STENCIL >
    >::type Add(const MULTI_INDEX_TYPE base_multi_index, const T_STENCIL_PROXY stencil_proxy);

    T& Center();
    const T& Center() const;
    T& operator()(const MULTI_INDEX_TYPE& multi_offset);
    const T& operator()(const MULTI_INDEX_TYPE& multi_offset) const;
    T& operator()(const int d, const int s);
    const T& operator()(const int d, const int s) const;

    template< class T_VECTOR >
    T Apply(const int base_index, const VECTOR<int,D>& strides, const T_VECTOR& x) const;
    template< class T_VECTOR >
    void Apply_Transpose(const int base_index, const VECTOR<int,D>& strides, T_VECTOR& y, const T x) const;
};

//#####################################################################
//#####################################################################

template< class T, int D >
inline CROSS_STENCIL<T,D>
CROSS_STENCIL<T,D>::
Construct_Zero()
{
    CROSS_STENCIL result;
    result.Zero();
    return result;
}

template< class T, int D >
inline void
CROSS_STENCIL<T,D>::
Zero()
{ std::fill(&values[0], &values[N_VALUES], static_cast<T>(0)); }

template< class T, int D >
inline int
CROSS_STENCIL<T,D>::
N_Zero() const
{ return static_cast< int >(std::count(&values[0], &values[N_VALUES], static_cast<T>(0))); }

template< class T, int D >
inline int
CROSS_STENCIL<T,D>::
N_Nonzero() const
{ return N_VALUES - N_Zero(); }

template< class T, int D >
inline T
CROSS_STENCIL<T,D>::
Sum() const
{ return std::accumulate(&values[0], &values[N_VALUES], static_cast<T>(0)); }

template< class T, int D >
inline CROSS_STENCIL<T,D>&
CROSS_STENCIL<T,D>::
operator+=(const CROSS_STENCIL& other)
{
    const T* p_srce = &other.values[0];
    T* p_dest = &values[0];
    PHYSBAM_LOOP( N_VALUES ) {
        *p_dest += *p_srce;
        ++p_srce;
        ++p_dest;
    }
    return *this;
}

template< class T, int D >
inline void
CROSS_STENCIL<T,D>::
Add(const VECTOR<T,D>& beta_dv_over_dx_dx, const MULTI_INDEX_TYPE& cell_multi_offset)
{
    for(int d = 1; d <= D; ++d) {
        assert(cell_multi_offset[d] == -1 || cell_multi_offset[d] == 0);
        const int s = 1 + 2 * cell_multi_offset[d];
        Center() += beta_dv_over_dx_dx[d];
        operator()(d,s) -= beta_dv_over_dx_dx[d];
    }
}

template< class T, int D >
template< class T_STENCIL_PROXY >
inline typename boost::enable_if<
    boost::is_same< typename T_STENCIL_PROXY::STENCIL_TYPE, CROSS_STENCIL<T,D> >
>::type
CROSS_STENCIL<T,D>::
Assign(const MULTI_INDEX_TYPE& /*base_multi_index*/, const T_STENCIL_PROXY& stencil_proxy)
{ operator=(stencil_proxy.Stencil()); }

template< class T, int D >
template< class T_STENCIL_PROXY >
inline typename boost::disable_if<
    boost::is_same< typename T_STENCIL_PROXY::STENCIL_TYPE, CROSS_STENCIL<T,D> >
>::type
CROSS_STENCIL<T,D>::
Assign(const MULTI_INDEX_TYPE base_multi_index, const T_STENCIL_PROXY stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, MULTI_INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    Zero();
    BOOST_FOREACH( typename T_STENCIL_PROXY::reference index_value, stencil_proxy )
        operator()(index_value.index - base_multi_index) = index_value.value;
}

template< class T, int D >
template< class T_STENCIL_PROXY >
inline typename boost::enable_if<
    boost::is_same< typename T_STENCIL_PROXY::STENCIL_TYPE, CROSS_STENCIL<T,D> >
>::type
CROSS_STENCIL<T,D>::
Add(const MULTI_INDEX_TYPE& /*base_multi_index*/, const T_STENCIL_PROXY& stencil_proxy)
{ operator+=(stencil_proxy.Stencil()); }

template< class T, int D >
template< class T_STENCIL_PROXY >
inline typename boost::disable_if<
    boost::is_same< typename T_STENCIL_PROXY::STENCIL_TYPE, CROSS_STENCIL<T,D> >
>::type
CROSS_STENCIL<T,D>::
Add(const MULTI_INDEX_TYPE base_multi_index, const T_STENCIL_PROXY stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, MULTI_INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    BOOST_FOREACH( typename T_STENCIL_PROXY::reference index_value, stencil_proxy )
        operator()(index_value.index - base_multi_index) += index_value.value;
}

template< class T, int D >
inline T&
CROSS_STENCIL<T,D>::
Center()
{ return values[D]; }

template< class T, int D >
inline const T&
CROSS_STENCIL<T,D>::
Center() const
{ return values[D]; }

template< class T, int D >
inline T&
CROSS_STENCIL<T,D>::
operator()(const MULTI_INDEX_TYPE& multi_offset)
{
    assert(multi_offset.L1_Norm() <= 1);
    for(int d = 1; d <= D; ++d)
        if(const int s = multi_offset[d])
            return operator()(d,s);
    return Center();
}

template< class T, int D >
inline const T&
CROSS_STENCIL<T,D>::
operator()(const MULTI_INDEX_TYPE& multi_offset) const
{
    assert(multi_offset.L1_Norm() <= 1);
    for(int d = 1; d <= D; ++d)
        if(const int s = multi_offset[d])
            return operator()(d,s);
    return Center();
}

template< class T, int D >
inline T&
CROSS_STENCIL<T,D>::
operator()(const int d, const int s)
{
    assert(1 <= d && d <= D);
    assert(s == -1 || s == +1);
    return values[D + s * (D - d + 1)];
}

template< class T, int D >
inline const T&
CROSS_STENCIL<T,D>::
operator()(const int d, const int s) const
{
    assert(1 <= d && d <= D);
    assert(s == -1 || s == +1);
    return values[D + s * (D - d + 1)];
}

template< class T, int D >
template< class T_VECTOR >
inline T
CROSS_STENCIL<T,D>::
Apply(const int base_index, const VECTOR<int,D>& strides, const T_VECTOR& x) const
{
    T y = Center() * x(base_index);
    for(int d = D; d >= 1; --d) {
        const int stride = strides[D - d + 1];
        if(values[D - d] != 0) {
            assert(1 <= base_index - stride);
            y += values[D - d] * x(base_index - stride);
        }
        if(values[D + d] != 0) {
            assert(base_index + stride <= x.Size());
            y += values[D + d] * x(base_index + stride);
        }
    }
    return y;
}

template< class T, int D >
template< class T_VECTOR >
inline void
CROSS_STENCIL<T,D>::
Apply_Transpose(const int base_index, const VECTOR<int,D>& strides, T_VECTOR& y, const T x) const
{
    y(base_index) = values[D] * x;
    for(int d = D; d >= 1; --d) {
        const int stride = strides[D - d + 1];
        if(values[D - d] != 0) {
            assert(1 <= base_index - stride);
            y(base_index - stride) += values[D - d] * x;
        }
        if(values[D + d] != 0) {
            assert(base_index + stride <= y.Size());
            y(base_index + stride) += values[D + d] * x;
        }
    }
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CROSS_STENCIL_HPP
