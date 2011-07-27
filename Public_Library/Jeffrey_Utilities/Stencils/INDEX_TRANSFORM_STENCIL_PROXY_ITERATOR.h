//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR_HPP

#include <boost/iterator/iterator_traits.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>

#include <Jeffrey_Utilities/ITERATOR_ADAPTOR.h>
#include <Jeffrey_Utilities/PROPAGATE_QUALIFIERS.h>
#include <Jeffrey_Utilities/RESULT_OF.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>

namespace PhysBAM
{

template< class T_BASE, class T_INDEX_TRANSFORM >
class INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR;

template< class T_BASE, class T_INDEX_TRANSFORM >
struct INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR_TRAITS
{
    typedef typename boost::iterator_value< T_BASE >::type BASE_INDEX_VALUE_TYPE;
    typedef typename BASE_INDEX_VALUE_TYPE::INDEX_TYPE BASE_INDEX_TYPE;
    typedef typename BASE_INDEX_VALUE_TYPE::VALUE_TYPE BASE_SCALAR_TYPE;
    typedef typename boost::iterator_reference< T_BASE >::type BASE_REFERENCE;
    typedef typename boost::remove_reference< BASE_REFERENCE >::type STRIPPED_BASE_REFERENCE_TYPE;
    typedef typename PROPAGATE_QUALIFIERS<
        BASE_REFERENCE,
        typename STRIPPED_BASE_REFERENCE_TYPE::INDEX_TYPE
    >::type BASE_REFERENCE_INDEX_TYPE;
    typedef typename PROPAGATE_QUALIFIERS<
        BASE_REFERENCE,
        typename STRIPPED_BASE_REFERENCE_TYPE::VALUE_TYPE
    >::type BASE_REFERENCE_SCALAR_TYPE;
    typedef typename RESULT_OF<
        const T_INDEX_TRANSFORM ( BASE_REFERENCE_INDEX_TYPE )
    >::type INDEX_TYPE;
    typedef BASE_SCALAR_TYPE SCALAR_TYPE;
    typedef BASE_REFERENCE_SCALAR_TYPE SCALAR_REFERENCE;
    typedef ITERATOR_ADAPTOR<
        INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE, T_INDEX_TRANSFORM >, // Derived
        T_BASE,                                      // Base
        INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE >,      // Value
        boost::use_default,                          // CategoryOrTraversal
        INDEX_VALUE< INDEX_TYPE, SCALAR_REFERENCE >, // Reference
        boost::use_default                           // Difference
    > ITERATOR_ADAPTOR_;
};

template< class T_BASE, class T_INDEX_TRANSFORM >
class INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR
    : public INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR_TRAITS< T_BASE, T_INDEX_TRANSFORM >::ITERATOR_ADAPTOR_
{
    typedef INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR_TRAITS< T_BASE, T_INDEX_TRANSFORM > INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR_TRAITS_;
    typedef typename INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR_TRAITS_::ITERATOR_ADAPTOR_ ITERATOR_ADAPTOR_;
public:
    typedef typename ITERATOR_ADAPTOR_::value_type value_type;
    typedef typename ITERATOR_ADAPTOR_::reference reference;
    typedef typename ITERATOR_ADAPTOR_::difference_type difference_type;

    T_INDEX_TRANSFORM index_transform;

    INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR();
    explicit INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR(const T_BASE& base_);
    INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR(
        const T_BASE& base_,
        const T_INDEX_TRANSFORM& index_transform_);

    template< class T_BASE2, class T_INDEX_TRANSFORM2 >
    INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR(
        const INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE2, T_INDEX_TRANSFORM2 >& other,
        typename boost::enable_if< boost::mpl::and_<
            boost::is_convertible< T_BASE2, T_BASE >,
            boost::is_convertible< T_INDEX_TRANSFORM2, T_INDEX_TRANSFORM >
        > >::type* = 0);

    using ITERATOR_ADAPTOR_::base;

    typedef typename INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR_TRAITS_::INDEX_TYPE INDEX_TYPE;
    typedef typename INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR_TRAITS_::SCALAR_TYPE SCALAR_TYPE;
    typedef typename INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR_TRAITS_::SCALAR_REFERENCE SCALAR_REFERENCE;
    INDEX_TYPE Index() const;
    SCALAR_REFERENCE Value() const;

private:
    friend class boost::iterator_core_access;

    reference dereference() const;
};

//#####################################################################
//#####################################################################

template< class T_BASE, class T_INDEX_TRANSFORM >
inline
INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE, T_INDEX_TRANSFORM >::
INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR()
{ }

template< class T_BASE, class T_INDEX_TRANSFORM >
inline
INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE, T_INDEX_TRANSFORM >::
INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR(const T_BASE& base_)
    : ITERATOR_ADAPTOR_(base_)
{ }

template< class T_BASE, class T_INDEX_TRANSFORM >
inline
INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE, T_INDEX_TRANSFORM >::
INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR(
    const T_BASE& base_,
    const T_INDEX_TRANSFORM& index_transform_)
    : ITERATOR_ADAPTOR_(base_),
      index_transform(index_transform_)
{ }

template< class T_BASE, class T_INDEX_TRANSFORM >
template< class T_BASE2, class T_INDEX_TRANSFORM2 >
inline
INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE, T_INDEX_TRANSFORM >::
INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR(
    const INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE2, T_INDEX_TRANSFORM2 >& other,
    typename boost::enable_if< boost::mpl::and_<
        boost::is_convertible< T_BASE2, T_BASE >,
        boost::is_convertible< T_INDEX_TRANSFORM2, T_INDEX_TRANSFORM >
    > >::type*)
    : ITERATOR_ADAPTOR_(other.base()),
      index_transform(other.index_transform)
{ }

template< class T_BASE, class T_INDEX_TRANSFORM >
inline typename INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE, T_INDEX_TRANSFORM >::INDEX_TYPE
INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE, T_INDEX_TRANSFORM >::
Index() const
{ return index_transform(base()->index); }

template< class T_BASE, class T_INDEX_TRANSFORM >
inline typename INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE, T_INDEX_TRANSFORM >::SCALAR_REFERENCE
INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE, T_INDEX_TRANSFORM >::
Value() const
{ return base()->value; }

template< class T_BASE, class T_INDEX_TRANSFORM >
inline typename INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE, T_INDEX_TRANSFORM >::reference
INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR< T_BASE, T_INDEX_TRANSFORM >::
dereference() const
{ return reference(Index(), Value()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_INDEX_TRANSFORM_STENCIL_PROXY_ITERATOR_HPP
