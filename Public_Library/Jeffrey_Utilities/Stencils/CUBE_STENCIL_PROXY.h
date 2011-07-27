//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CUBE_STENCIL_PROXY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CUBE_STENCIL_PROXY_HPP

#include <boost/concept/assert.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/utility/enable_if.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/PROPAGATE_CONST.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_ITERATOR.h>

namespace PhysBAM
{

template< class T_STENCIL >
struct CUBE_STENCIL_PROXY
{
    typedef typename T_STENCIL::MULTI_INDEX_TYPE INDEX_TYPE;
    typedef typename T_STENCIL::SCALAR_TYPE SCALAR_TYPE;
    typedef INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE > INDEX_VALUE_TYPE;

    static const int DIMENSION = T_STENCIL::DIMENSION;
    static const int MIN_OFFSET = T_STENCIL::MIN_OFFSET;
    static const int MAX_OFFSET = T_STENCIL::MAX_OFFSET;
    static const int WIDTH = T_STENCIL::WIDTH;

    typedef typename T_STENCIL::MULTI_INDEX_TYPE MULTI_INDEX_TYPE;
    typedef MULTI_INDEX_CUBE< DIMENSION, MIN_OFFSET, MAX_OFFSET > MULTI_INDEX_CUBE_TYPE;

    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS(
        CUBE_STENCIL_PROXY,
        (( typename MULTI_INDEX_CUBE_TYPE const, multi_index_cube ))
        (( typename T_STENCIL&, stencil ))
    )

    CUBE_STENCIL_PROXY(const MULTI_INDEX_TYPE& base_multi_index, T_STENCIL& stencil_);

    template< class T_STENCIL2 >
    CUBE_STENCIL_PROXY(const CUBE_STENCIL_PROXY< T_STENCIL2 >& other,
        typename boost::enable_if< boost::is_convertible< T_STENCIL2*, T_STENCIL* > >::type* = 0);

    CUBE_STENCIL_PROXY& operator+=(const INDEX_VALUE_TYPE& index_value);

    CUBE_STENCIL_PROXY& operator=(const CUBE_STENCIL_PROXY& other);
    template< class T_STENCIL_PROXY >
    CUBE_STENCIL_PROXY& operator=(const T_STENCIL_PROXY& stencil_proxy);
    template< class T_STENCIL_PROXY >
    CUBE_STENCIL_PROXY& operator+=(const T_STENCIL_PROXY& stencil_proxy);

    int N_Nonzero() const;

    typedef STENCIL_PROXY_ITERATOR<
        typename MULTI_INDEX_CUBE_TYPE::iterator,
        typename PROPAGATE_CONST< T_STENCIL, SCALAR_TYPE >::type *
    > iterator;
    typedef iterator const_iterator;
    typedef typename iterator::reference reference;
    iterator begin() const;
    iterator end() const;
};

//#####################################################################
//#####################################################################

template< class T_STENCIL >
inline
CUBE_STENCIL_PROXY< T_STENCIL >::
CUBE_STENCIL_PROXY(const MULTI_INDEX_TYPE& base_multi_index, T_STENCIL& stencil_)
    : multi_index_cube(base_multi_index),
      stencil(stencil_)
{ }

template< class T_STENCIL >
template< class T_STENCIL2 >
CUBE_STENCIL_PROXY< T_STENCIL >::
CUBE_STENCIL_PROXY(const CUBE_STENCIL_PROXY< T_STENCIL2 >& other,
    typename boost::enable_if< boost::is_convertible< T_STENCIL2*, T_STENCIL* > >::type*)
    : multi_index_cube(other.multi_index_cube),
      stencil(other.stencil)
{ }

template< class T_STENCIL >
inline CUBE_STENCIL_PROXY< T_STENCIL >&
CUBE_STENCIL_PROXY< T_STENCIL >::
operator+=(const INDEX_VALUE_TYPE& index_value)
{
    stencil(multi_index_cube.Multi_Offset(index_value.index)) += index_value.value;
    return *this;
}

template< class T_STENCIL >
inline CUBE_STENCIL_PROXY< T_STENCIL >&
CUBE_STENCIL_PROXY< T_STENCIL >::
operator=(const CUBE_STENCIL_PROXY& other)
{ return operator=< CUBE_STENCIL_PROXY >(other); }

template< class T_STENCIL >
template< class T_STENCIL_PROXY >
inline CUBE_STENCIL_PROXY< T_STENCIL >&
CUBE_STENCIL_PROXY< T_STENCIL >::
operator=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    stencil.Assign(multi_index_cube.base_multi_index, stencil_proxy);
    return *this;
}

template< class T_STENCIL >
template< class T_STENCIL_PROXY >
inline CUBE_STENCIL_PROXY< T_STENCIL >&
CUBE_STENCIL_PROXY< T_STENCIL >::
operator+=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    stencil.Add(multi_index_cube.base_multi_index, stencil_proxy);
    return *this;
}

template< class T_STENCIL >
inline int
CUBE_STENCIL_PROXY< T_STENCIL >::
N_Nonzero() const
{ return stencil.N_Nonzero(); }

template< class T_STENCIL >
inline typename CUBE_STENCIL_PROXY< T_STENCIL >::iterator
CUBE_STENCIL_PROXY< T_STENCIL >::
begin() const
{ return iterator(multi_index_cube.begin(), &stencil.values[0]); }

template< class T_STENCIL >
inline typename CUBE_STENCIL_PROXY< T_STENCIL >::iterator
CUBE_STENCIL_PROXY< T_STENCIL >::
end() const
{ return iterator(multi_index_cube.end(), &stencil.values[T_STENCIL::N_VALUES]); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CUBE_STENCIL_PROXY_HPP
