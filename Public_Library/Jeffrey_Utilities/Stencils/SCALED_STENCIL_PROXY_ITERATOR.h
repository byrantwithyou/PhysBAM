//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SCALED_STENCIL_PROXY_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SCALED_STENCIL_PROXY_ITERATOR_HPP

#include <boost/iterator/iterator_traits.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/utility/enable_if.hpp>

#include <Jeffrey_Utilities/ITERATOR_ADAPTOR.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>

namespace PhysBAM
{

template< class T_BASE >
class SCALED_STENCIL_PROXY_ITERATOR;

template< class T_BASE >
struct SCALED_STENCIL_PROXY_ITERATOR_TRAITS
{
    typedef typename boost::iterator_value< T_BASE >::type::INDEX_TYPE INDEX_TYPE;
    typedef typename boost::iterator_value< T_BASE >::type::VALUE_TYPE SCALAR_TYPE;
    typedef typename boost::remove_reference<
        typename boost::iterator_reference< T_BASE >::type
    >::type::INDEX_TYPE INDEX_REFERENCE_TYPE;
    typedef ITERATOR_ADAPTOR<
        SCALED_STENCIL_PROXY_ITERATOR< T_BASE >,          // Derived
        T_BASE,                                           // Base
        INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE >,           // Value
        boost::use_default,                               // CategoryOrTraversal
        INDEX_VALUE< INDEX_REFERENCE_TYPE, SCALAR_TYPE >, // Reference
        boost::use_default                                // Difference
    > ITERATOR_ADAPTOR_;
};

template< class T_BASE >
class SCALED_STENCIL_PROXY_ITERATOR
    : public SCALED_STENCIL_PROXY_ITERATOR_TRAITS< T_BASE >::ITERATOR_ADAPTOR_
{
    typedef SCALED_STENCIL_PROXY_ITERATOR_TRAITS< T_BASE > SCALED_STENCIL_PROXY_ITERATOR_TRAITS_;
    typedef typename SCALED_STENCIL_PROXY_ITERATOR_TRAITS_::ITERATOR_ADAPTOR_ ITERATOR_ADAPTOR_;
public:
    typedef typename ITERATOR_ADAPTOR_::value_type value_type;
    typedef typename ITERATOR_ADAPTOR_::reference reference;
    typedef typename ITERATOR_ADAPTOR_::difference_type difference_type;

    typedef typename SCALED_STENCIL_PROXY_ITERATOR_TRAITS_::INDEX_TYPE INDEX_TYPE;
    typedef typename SCALED_STENCIL_PROXY_ITERATOR_TRAITS_::SCALAR_TYPE SCALAR_TYPE;

    SCALED_STENCIL_PROXY_ITERATOR();
    explicit SCALED_STENCIL_PROXY_ITERATOR(const T_BASE& base);
    SCALED_STENCIL_PROXY_ITERATOR(const T_BASE& base, const SCALAR_TYPE scaling);

    template< class T_BASE2 >
    SCALED_STENCIL_PROXY_ITERATOR(const SCALED_STENCIL_PROXY_ITERATOR< T_BASE2 >& other,
        typename boost::enable_if< boost::is_convertible< T_BASE2, T_BASE > >::type* = 0);

    using ITERATOR_ADAPTOR_::base;

    INDEX_TYPE Index() const;
    SCALAR_TYPE Value() const;
    SCALAR_TYPE Scaling() const;

private:
    SCALAR_TYPE m_scaling;

    friend class boost::iterator_core_access;

    reference dereference() const;
};

//#####################################################################
//#####################################################################

template< class T_BASE >
inline
SCALED_STENCIL_PROXY_ITERATOR< T_BASE >::
SCALED_STENCIL_PROXY_ITERATOR()
{ }

template< class T_BASE >
inline
SCALED_STENCIL_PROXY_ITERATOR< T_BASE >::
SCALED_STENCIL_PROXY_ITERATOR(const T_BASE& base)
    : ITERATOR_ADAPTOR_(base),
      m_scaling(1)
{ }

template< class T_BASE >
inline
SCALED_STENCIL_PROXY_ITERATOR< T_BASE >::
SCALED_STENCIL_PROXY_ITERATOR(const T_BASE& base, const SCALAR_TYPE scaling)
    : ITERATOR_ADAPTOR_(base),
      m_scaling(scaling)
{ }

template< class T_BASE >
template< class T_BASE2 >
inline
SCALED_STENCIL_PROXY_ITERATOR< T_BASE >::
SCALED_STENCIL_PROXY_ITERATOR(const SCALED_STENCIL_PROXY_ITERATOR< T_BASE2 >& other,
    typename boost::enable_if< boost::is_convertible< T_BASE2, T_BASE > >::type*)
    : ITERATOR_ADAPTOR_(other.base()),
      m_scaling(other.Scaling())
{ }

template< class T_BASE >
inline typename SCALED_STENCIL_PROXY_ITERATOR< T_BASE >::INDEX_TYPE
SCALED_STENCIL_PROXY_ITERATOR< T_BASE >::
Index() const
{ return base()->index; }

template< class T_BASE >
inline typename SCALED_STENCIL_PROXY_ITERATOR< T_BASE >::SCALAR_TYPE
SCALED_STENCIL_PROXY_ITERATOR< T_BASE >::
Value() const
{ return m_scaling * base()->value; }

template< class T_BASE >
inline typename SCALED_STENCIL_PROXY_ITERATOR< T_BASE >::SCALAR_TYPE
SCALED_STENCIL_PROXY_ITERATOR< T_BASE >::
Scaling() const
{ return m_scaling; }

template< class T_BASE >
inline typename SCALED_STENCIL_PROXY_ITERATOR< T_BASE >::reference
SCALED_STENCIL_PROXY_ITERATOR< T_BASE >::
dereference() const
{ return m_scaling * *base(); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_SCALED_STENCIL_PROXY_ITERATOR_HPP
