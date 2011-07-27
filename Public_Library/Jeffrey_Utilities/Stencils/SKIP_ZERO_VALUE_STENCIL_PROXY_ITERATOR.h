//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR_HPP

#include <cassert>

#include <boost/iterator/iterator_categories.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>

#include <Jeffrey_Utilities/ITERATOR_ADAPTOR.h>

namespace PhysBAM
{

template< class T_BASE >
class SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR;

template< class T_BASE >
struct SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR_TRAITS
{
    typedef ITERATOR_ADAPTOR<
        SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR< T_BASE >, // Derived
        T_BASE,                       // Base
        boost::use_default,           // Value
        boost::forward_traversal_tag, // CategoryOrTraversal
        boost::use_default,           // Reference
        boost::use_default            // Difference
    > ITERATOR_ADAPTOR_;
};

template< class T_BASE >
class SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR
    : public SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR_TRAITS< T_BASE >::ITERATOR_ADAPTOR_
{
    typedef SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR_TRAITS< T_BASE > SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR_TRAITS_;
    typedef typename SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR_TRAITS_::ITERATOR_ADAPTOR_ ITERATOR_ADAPTOR_;
public:
    typedef typename ITERATOR_ADAPTOR_::value_type value_type;
    typedef typename ITERATOR_ADAPTOR_::reference reference;
    typedef typename ITERATOR_ADAPTOR_::difference_type difference_type;

    SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR();
    SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR(
        const T_BASE& base_,
        const T_BASE& base_end);

    template< class T_BASE2 >
    SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR(
        const SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR< T_BASE2 >& other,
        typename boost::enable_if<
            boost::is_convertible< T_BASE2, T_BASE >
        >::type* = 0);

    using ITERATOR_ADAPTOR_::base;
    const T_BASE& Base_End() const;

private:
    T_BASE m_base_end;

    using ITERATOR_ADAPTOR_::base_reference;

    friend class boost::iterator_core_access;
    void increment();
};

//#####################################################################
//#####################################################################

template< class T_BASE >
inline
SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR< T_BASE >::
SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR()
{ }

template< class T_BASE >
inline
SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR< T_BASE >::
SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR(
    const T_BASE& base_,
    const T_BASE& base_end)
    : ITERATOR_ADAPTOR_(base_),
      m_base_end(base_end)
{ for(; base() != m_base_end && base()->value == 0; ++base_reference()); }

template< class T_BASE >
template< class T_BASE2 >
inline
SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR< T_BASE >::
SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR(
    const SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR< T_BASE2 >& other,
    typename boost::enable_if<
        boost::is_convertible< T_BASE2, T_BASE >
    >::type*)
    : ITERATOR_ADAPTOR_(other.base()),
      m_base_end(other.Base_End())
{ }

template< class T_BASE >
inline void
SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR< T_BASE >::
increment()
{
    assert(base() != m_base_end);
    while(++base_reference() != m_base_end && base()->value == 0);
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SKIP_ZERO_VALUE_STENCIL_PROXY_ITERATOR_HPP
