//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CUBE_STENCIL_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CUBE_STENCIL_HPP

#include <cassert>

#include <algorithm>
#include <numeric>

#include <boost/concept/assert.hpp>
#include <boost/foreach.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/LOOP.h>
#include <Jeffrey_Utilities/Math/STATIC_MIN.h>
#include <Jeffrey_Utilities/Math/STATIC_POW.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_FWD.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>

namespace PhysBAM
{

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ /* = -MIN_OFFSET_ */ >
struct CUBE_STENCIL
{
    BOOST_MPL_ASSERT_RELATION( MIN_OFFSET_, <=, MAX_OFFSET_ );

    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef T SCALAR_TYPE;

    static const int DIMENSION = D;
    static const int MIN_OFFSET = MIN_OFFSET_;
    static const int MAX_OFFSET = MAX_OFFSET_;
    static const int WIDTH = 1 + MAX_OFFSET - MIN_OFFSET;
    static const int N_VALUES = STATIC_POW_C< WIDTH, D >::value;

    T values[N_VALUES];

    static CUBE_STENCIL Construct_Zero();

    void Zero();
    int N_Zero() const;
    int N_Nonzero() const;

    template< class T_STENCIL_PROXY >
    void Assign(
        const MULTI_INDEX_TYPE base_multi_index,
        const T_STENCIL_PROXY stencil_proxy);

    template< class T_STENCIL_PROXY >
    void Add(
        const MULTI_INDEX_TYPE base_multi_index,
        const T_STENCIL_PROXY stencil_proxy);

    T& operator()(const MULTI_INDEX_TYPE& multi_offset);
    const T& operator()(const MULTI_INDEX_TYPE& multi_offset) const;

    T Sum() const;

    template< class T_VECTOR >
    T Apply(const int base_index, const VECTOR<int,D>& strides, const T_VECTOR& x) const;
    template< class T_VECTOR >
    void Apply_Transpose(const int base_index, const VECTOR<int,D>& strides, T_VECTOR& y, const T x) const;
};

//#####################################################################
//#####################################################################

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline CUBE_STENCIL< T, D, MIN_OFFSET_, MAX_OFFSET_ >
CUBE_STENCIL< T, D, MIN_OFFSET_, MAX_OFFSET_ >::
Construct_Zero()
{
    CUBE_STENCIL result;
    result.Zero();
    return result;
}

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline void
CUBE_STENCIL< T, D, MIN_OFFSET_, MAX_OFFSET_ >::
Zero()
{ std::fill(&values[0], &values[N_VALUES], static_cast<T>(0)); }

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline int
CUBE_STENCIL< T, D, MIN_OFFSET_, MAX_OFFSET_ >::
N_Zero() const
{ return static_cast< int >(std::count(&values[0], &values[N_VALUES], static_cast<T>(0))); }

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline int
CUBE_STENCIL< T, D, MIN_OFFSET_, MAX_OFFSET_ >::
N_Nonzero() const
{ return N_VALUES - N_Zero(); }

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ >
template< class T_STENCIL_PROXY >
inline void
CUBE_STENCIL< T, D, MIN_OFFSET_, MAX_OFFSET_ >::
Assign(
    const MULTI_INDEX_TYPE base_multi_index,
    const T_STENCIL_PROXY stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, MULTI_INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    Zero();
    BOOST_FOREACH( typename T_STENCIL_PROXY::reference index_value, stencil_proxy )
        operator()(index_value.index - base_multi_index) = index_value.value;
}

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ >
template< class T_STENCIL_PROXY >
inline void
CUBE_STENCIL< T, D, MIN_OFFSET_, MAX_OFFSET_ >::
Add(
    const MULTI_INDEX_TYPE base_multi_index,
    const T_STENCIL_PROXY stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, MULTI_INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    BOOST_FOREACH( typename T_STENCIL_PROXY::reference index_value, stencil_proxy )
        operator()(index_value.index - base_multi_index) += index_value.value;
}

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline T&
CUBE_STENCIL< T, D, MIN_OFFSET_, MAX_OFFSET_ >::
operator()(const MULTI_INDEX_TYPE& multi_offset)
{ return values[STATIC_MULTI_INDEX_CUBE< D, MIN_OFFSET, MAX_OFFSET >::Linear_Index(multi_offset) - 1]; }

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline const T&
CUBE_STENCIL< T, D, MIN_OFFSET_, MAX_OFFSET_ >::
operator()(const MULTI_INDEX_TYPE& multi_offset) const
{ return values[STATIC_MULTI_INDEX_CUBE< D, MIN_OFFSET, MAX_OFFSET >::Linear_Index(multi_offset) - 1]; }

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline T
CUBE_STENCIL< T, D, MIN_OFFSET_, MAX_OFFSET_ >::
Sum() const
{ return std::accumulate(&values[0], &values[N_VALUES], static_cast<T>(0)); }

namespace Detail_CUBE_STENCIL
{

template< class T, int D, int MIN_OFFSET, int MAX_OFFSET >
struct DISPATCH;

template< class T, int MIN_OFFSET, int MAX_OFFSET >
struct DISPATCH< T, 1, MIN_OFFSET, MAX_OFFSET >
{
    static const int WIDTH = 1 + MAX_OFFSET - MIN_OFFSET;
    template< class T_VECTOR >
    static T Apply(const T* p_value, int index, VECTOR<int,1> strides, const T_VECTOR& x)
    {
        index += MIN_OFFSET * strides.Sum();
        T y = 0;
        PHYSBAM_LOOP( WIDTH ) {
            if(*p_value != 0) {
                assert(1 <= index && index <= x.Size());
                y += *p_value * x(index);
            }
            ++p_value;
            index += strides[1];
        }
        return y;
    }
    template< class T_VECTOR >
    static void Apply_Transpose(const T* p_value, int index, VECTOR<int,1> strides, T_VECTOR& y, const T x)
    {
        index += MIN_OFFSET * strides.Sum();
        PHYSBAM_LOOP( WIDTH ) {
            if(*p_value != 0) {
                assert(1 <= index && index <= y.Size());
                y(index) += *p_value * x;
            }
            ++p_value;
            index += strides[1];
        }
    }
};

template< class T, int MIN_OFFSET, int MAX_OFFSET >
struct DISPATCH< T, 2, MIN_OFFSET, MAX_OFFSET >
{
    static const int WIDTH = 1 + MAX_OFFSET - MIN_OFFSET;
    template< class T_VECTOR >
    static T Apply(const T* p_value, int index, VECTOR<int,2> strides, const T_VECTOR& x)
    {
        index += MIN_OFFSET * strides.Sum();
        strides[1] -= WIDTH * strides[2];
        T y = 0;
        PHYSBAM_LOOP( WIDTH ) {
            PHYSBAM_LOOP( WIDTH ) {
                if(*p_value != 0) {
                    assert(1 <= index && index <= x.Size());
                    y += *p_value * x(index);
                }
                ++p_value;
                index += strides[2];
            }
            index += strides[1];
        }
        return y;
    }
    template< class T_VECTOR >
    static void Apply_Transpose(const T* p_value, int index, VECTOR<int,2> strides, T_VECTOR& y, const T x)
    {
        index += MIN_OFFSET * strides.Sum();
        strides[1] -= WIDTH * strides[2];
        PHYSBAM_LOOP( WIDTH ) {
            PHYSBAM_LOOP( WIDTH ) {
                if(*p_value != 0) {
                    assert(1 <= index && index <= y.Size());
                    y(index) += *p_value * x;
                }
                ++p_value;
                index += strides[2];
            }
            index += strides[1];
        }
    }
};

template< class T, int MIN_OFFSET, int MAX_OFFSET >
struct DISPATCH< T, 3, MIN_OFFSET, MAX_OFFSET >
{
    static const int WIDTH = 1 + MAX_OFFSET - MIN_OFFSET;
    template< class T_VECTOR >
    static T Apply(const T* p_value, int index, VECTOR<int,3> strides, const T_VECTOR& x)
    {
        index += MIN_OFFSET * strides.Sum();
        strides[1] -= WIDTH * strides[2];
        strides[2] -= WIDTH * strides[3];
        T y = 0;
        PHYSBAM_LOOP( WIDTH ) {
            PHYSBAM_LOOP( WIDTH ) {
                PHYSBAM_LOOP( WIDTH ) {
                    if(*p_value != 0) {
                        assert(1 <= index && index <= x.Size());
                        y += *p_value * x(index);
                    }
                    ++p_value;
                    index += strides[3];
                }
                index += strides[2];
            }
            index += strides[1];
        }
        return y;
    }
    template< class T_VECTOR >
    static void Apply_Transpose(const T* p_value, int index, VECTOR<int,3> strides, T_VECTOR& y, const T x)
    {
        index += MIN_OFFSET * strides.Sum();
        strides[1] -= WIDTH * strides[2];
        strides[2] -= WIDTH * strides[3];
        PHYSBAM_LOOP( WIDTH ) {
            PHYSBAM_LOOP( WIDTH ) {
                PHYSBAM_LOOP( WIDTH ) {
                    if(*p_value != 0) {
                        assert(1 <= index && index <= y.Size());
                        y(index) += *p_value * x;
                    }
                    ++p_value;
                    index += strides[3];
                }
                index += strides[2];
            }
            index += strides[1];
        }
    }
};

} // namespace Detail_CUBE_STENCIL

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ >
template< class T_VECTOR >
inline T
CUBE_STENCIL< T, D, MIN_OFFSET_, MAX_OFFSET_ >::
Apply(const int base_index, const VECTOR<int,D>& strides, const T_VECTOR& x) const
{ return Detail_CUBE_STENCIL::DISPATCH< T, D, MIN_OFFSET, MAX_OFFSET >::Apply(&values[0], base_index, strides, x); }

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ >
template< class T_VECTOR >
inline void
CUBE_STENCIL< T, D, MIN_OFFSET_, MAX_OFFSET_ >::
Apply_Transpose(const int base_index, const VECTOR<int,D>& strides, T_VECTOR& y, const T x) const
{ Detail_CUBE_STENCIL::DISPATCH< T, D, MIN_OFFSET, MAX_OFFSET >::Apply_Transpose(&values[0], base_index, strides, y, x); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CUBE_STENCIL_HPP
