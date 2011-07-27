//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_INDEX_TRANSFORM_STENCIL_PROXY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_INDEX_TRANSFORM_STENCIL_PROXY_HPP

#include <boost/concept/assert.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/REMOVE_QUALIFIERS.h>
#include <Jeffrey_Utilities/RESULT_OF.h>
#include <Jeffrey_Utilities/Stencils/INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_CONCEPT.h>
#include <Jeffrey_Utilities/Stencils/STENCIL_PROXY_FWD.h>

namespace PhysBAM
{

template<
    class T_BASE,
    class T_INDEX_TRANSFORM,
    class T_INVERSE_INDEX_TRANSFORM /*= T_INDEX_TRANSFORM*/
>
struct INDEX_TRANSFORM_STENCIL_PROXY
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_BASE >));

    typedef INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR<
        typename T_BASE::iterator, T_INDEX_TRANSFORM
    > iterator;

    typedef typename iterator::INDEX_TYPE INDEX_TYPE;
    typedef typename iterator::SCALAR_TYPE SCALAR_TYPE;
    typedef INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE > INDEX_VALUE_TYPE;

    BOOST_MPL_ASSERT((boost::is_convertible<
        typename RESULT_OF<
            const T_INVERSE_INDEX_TRANSFORM (
                typename RESULT_OF<
                    const T_INDEX_TRANSFORM (
                        typename T_BASE::INDEX_TYPE
                    )
                >::type
            )
        >::type,
        typename T_BASE::INDEX_TYPE
    >));
    BOOST_MPL_ASSERT((boost::is_convertible<
        typename RESULT_OF<
            const T_INDEX_TRANSFORM (
                typename RESULT_OF<
                    const T_INVERSE_INDEX_TRANSFORM ( INDEX_TYPE )
                >::type
            )
        >::type,
        INDEX_TYPE
    >));

    typedef T_BASE BASE_TYPE;
    typedef T_INDEX_TRANSFORM INDEX_TRANSFORM_TYPE;
    typedef T_INVERSE_INDEX_TRANSFORM INVERSE_INDEX_TRANSFORM_TYPE;

    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS(
        INDEX_TRANSFORM_STENCIL_PROXY,
        (( typename T_BASE, base ))
        (( typename T_INDEX_TRANSFORM const, index_transform ))
        (( typename T_INVERSE_INDEX_TRANSFORM const, inverse_index_transform ))
    )

    INDEX_TRANSFORM_STENCIL_PROXY(
        const T_BASE& base_,
        const T_INDEX_TRANSFORM& index_transform);

    template< class T_BASE2, class T_INDEX_TRANSFORM2, class T_INVERSE_INDEX_TRANSFORM2 >
    INDEX_TRANSFORM_STENCIL_PROXY(
        const INDEX_TRANSFORM_STENCIL_PROXY< T_BASE2, T_INDEX_TRANSFORM2, T_INVERSE_INDEX_TRANSFORM2 >& other,
        typename boost::enable_if< boost::mpl::and_<
            boost::is_convertible< T_BASE2, T_BASE >,
            boost::is_convertible< T_INDEX_TRANSFORM2, T_INDEX_TRANSFORM >,
            boost::is_convertible< T_INVERSE_INDEX_TRANSFORM2, T_INVERSE_INDEX_TRANSFORM >
        > >::type* = 0);

    INDEX_TRANSFORM_STENCIL_PROXY& operator+=(const INDEX_VALUE_TYPE& index_value);

    template< class T_STENCIL_PROXY >
    INDEX_TRANSFORM_STENCIL_PROXY& operator=(const T_STENCIL_PROXY& stencil_proxy);
    template< class T_STENCIL_PROXY >
    INDEX_TRANSFORM_STENCIL_PROXY& operator+=(const T_STENCIL_PROXY& stencil_proxy);

    int N_Nonzero() const;

    typedef iterator const_iterator;
    typedef typename iterator::reference reference;
    iterator begin() const;
    iterator end() const;
};

template< class T_BASE, class T_INDEX_TRANSFORM, class T_INVERSE_INDEX_TRANSFORM >
inline INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >
Make_Index_Transform_Stencil_Proxy(
    const T_BASE& base,
    const T_INDEX_TRANSFORM& index_transform,
    const T_INVERSE_INDEX_TRANSFORM& inverse_index_transform)
{
    return INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >(
        base, index_transform, inverse_index_transform
    );
}

template< class T_BASE, class T_INDEX_TRANSFORM >
inline INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM >
Make_Index_Transform_Stencil_Proxy(
    const T_BASE& base,
    const T_INDEX_TRANSFORM& index_transform)
{
    return INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM >(
        base, index_transform
    );
}

template< class T_INDEX_TRANSFORM, class T_INVERSE_INDEX_TRANSFORM = T_INDEX_TRANSFORM >
struct INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION,
        (( typename T_INDEX_TRANSFORM const, index_transform ))
        (( typename T_INVERSE_INDEX_TRANSFORM const, inverse_index_transform ))
    )
public:
    explicit INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION(const T_INDEX_TRANSFORM& index_transform_)
        : index_transform(index_transform_),
          inverse_index_transform(index_transform_)
    { }

    template<class> struct result;
    template< class T_THIS, class T_BASE >
    struct result< T_THIS ( T_BASE ) >
    {
        typedef INDEX_TRANSFORM_STENCIL_PROXY<
            typename REMOVE_QUALIFIERS< T_BASE >::type,
            T_INDEX_TRANSFORM,
            T_INVERSE_INDEX_TRANSFORM
        > type;
    };

    template< class T_BASE >
    INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >
    operator()(const T_BASE& base) const
    {
        return INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >(
            base, index_transform, inverse_index_transform
        );
    }
};

template< class T_INDEX_TRANSFORM, class T_INVERSE_INDEX_TRANSFORM >
inline INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION< T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >
Make_Index_Transform_Stencil_Proxy_Function(
    const T_INDEX_TRANSFORM& index_transform,
    const T_INVERSE_INDEX_TRANSFORM& inverse_index_transform)
{
    return INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION< T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >(
        index_transform, inverse_index_transform
    );
}

template< class T_INDEX_TRANSFORM >
inline INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION< T_INDEX_TRANSFORM >
Make_Index_Transform_Stencil_Proxy_Function(const T_INDEX_TRANSFORM& index_transform)
{ return INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION< T_INDEX_TRANSFORM >(index_transform); }

//#####################################################################
//#####################################################################

template< class T_BASE, class T_INDEX_TRANSFORM, class T_INVERSE_INDEX_TRANSFORM >
inline
INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >::
INDEX_TRANSFORM_STENCIL_PROXY(
    const T_BASE& base_,
    const T_INDEX_TRANSFORM& index_transform)
    : base(base_),
      index_transform(index_transform),
      inverse_index_transform(index_transform)
{ }

template< class T_BASE, class T_INDEX_TRANSFORM, class T_INVERSE_INDEX_TRANSFORM >
template< class T_BASE2, class T_INDEX_TRANSFORM2, class T_INVERSE_INDEX_TRANSFORM2 >
inline
INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >::
INDEX_TRANSFORM_STENCIL_PROXY(
    const INDEX_TRANSFORM_STENCIL_PROXY< T_BASE2, T_INDEX_TRANSFORM2, T_INVERSE_INDEX_TRANSFORM2 >& other,
    typename boost::enable_if< boost::mpl::and_<
        boost::is_convertible< T_BASE2, T_BASE >,
        boost::is_convertible< T_INDEX_TRANSFORM2, T_INDEX_TRANSFORM >,
        boost::is_convertible< T_INVERSE_INDEX_TRANSFORM2, T_INVERSE_INDEX_TRANSFORM >
    > >::type*)
    : base(other.base),
      index_transform(other.index_transform),
      inverse_index_transform(other.inverse_index_transform)
{ }

template< class T_BASE, class T_INDEX_TRANSFORM, class T_INVERSE_INDEX_TRANSFORM >
inline INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >&
INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >::
operator+=(const INDEX_VALUE_TYPE& index_value)
{
    typedef INDEX_VALUE<
        typename T_BASE::INDEX_TYPE,
        SCALAR_TYPE
    > BASE_INDEX_VALUE_TYPE;
    base += BASE_INDEX_VALUE_TYPE(
        index_transform(index_value.index),
        index_value.value
    );
    return *this;
}

template< class T_BASE, class T_INDEX_TRANSFORM, class T_INVERSE_INDEX_TRANSFORM >
template< class T_STENCIL_PROXY >
inline INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >&
INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >::
operator=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    base = Make_Index_Transform_Stencil_Proxy(stencil_proxy, inverse_index_transform, index_transform);
    return *this;
}

template< class T_BASE, class T_INDEX_TRANSFORM, class T_INVERSE_INDEX_TRANSFORM >
template< class T_STENCIL_PROXY >
inline INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >&
INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >::
operator+=(const T_STENCIL_PROXY& stencil_proxy)
{
    BOOST_CONCEPT_ASSERT((CONST_STENCIL_PROXY_CONCEPT< T_STENCIL_PROXY >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::INDEX_TYPE, INDEX_TYPE >));
    BOOST_MPL_ASSERT((boost::is_same< typename T_STENCIL_PROXY::SCALAR_TYPE, SCALAR_TYPE >));
    base += Make_Index_Transform_Stencil_Proxy(stencil_proxy, inverse_index_transform, index_transform);
    return *this;
}

template< class T_BASE, class T_INDEX_TRANSFORM, class T_INVERSE_INDEX_TRANSFORM >
inline int
INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >::
N_Nonzero() const
{ return base.N_Nonzero(); }

template< class T_BASE, class T_INDEX_TRANSFORM, class T_INVERSE_INDEX_TRANSFORM >
inline typename INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >::iterator
INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >::
begin() const
{ return iterator(base.begin(), index_transform); }

template< class T_BASE, class T_INDEX_TRANSFORM, class T_INVERSE_INDEX_TRANSFORM >
inline typename INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >::iterator
INDEX_TRANSFORM_STENCIL_PROXY< T_BASE, T_INDEX_TRANSFORM, T_INVERSE_INDEX_TRANSFORM >::
end() const
{ return iterator(base.end(), index_transform); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_INDEX_TRANSFORM_STENCIL_PROXY_HPP
