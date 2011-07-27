//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_STENCIL_PROXY_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_STENCIL_PROXY_ITERATOR_HPP

#include <cassert>

#include <iterator>

#include <boost/iterator/iterator_traits.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/not.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_reference.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/utility/enable_if.hpp>

#include <Jeffrey_Utilities/ITERATOR_FACADE.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>

namespace PhysBAM
{

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
class STENCIL_PROXY_ITERATOR;

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
struct STENCIL_PROXY_ITERATOR_TRAITS
{
    typedef typename boost::iterator_value< T_INDEX_ITERATOR >::type INDEX_TYPE;
    typedef typename boost::iterator_reference< T_INDEX_ITERATOR >::type INDEX_REFERENCE_TYPE;
    typedef typename boost::iterator_value< T_VALUE_ITERATOR >::type VALUE_TYPE;
    typedef typename boost::iterator_reference< T_VALUE_ITERATOR >::type VALUE_REFERENCE_TYPE;
    typedef INDEX_VALUE< INDEX_TYPE, VALUE_TYPE > INDEX_VALUE_TYPE;
    typedef typename boost::mpl::if_<
        boost::mpl::and_<
            boost::is_reference< VALUE_REFERENCE_TYPE >,
            boost::mpl::not_< boost::is_const<
                typename boost::remove_reference< VALUE_REFERENCE_TYPE >::type
            > >
        >,
        INDEX_VALUE_TYPE,
        const INDEX_VALUE_TYPE
    >::type FACADE_INDEX_VALUE_TYPE;
    typedef ITERATOR_FACADE<
        STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >, // Derived
        FACADE_INDEX_VALUE_TYPE,                                      // Value
        std::random_access_iterator_tag,                              // CategoryOrTraversal
        INDEX_VALUE< INDEX_TYPE, VALUE_REFERENCE_TYPE >,              // Reference
        int                                                           // Difference
    > ITERATOR_FACADE_;
};

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
class STENCIL_PROXY_ITERATOR
    : public STENCIL_PROXY_ITERATOR_TRAITS< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::ITERATOR_FACADE_
{
    typedef STENCIL_PROXY_ITERATOR_TRAITS< T_INDEX_ITERATOR, T_VALUE_ITERATOR > STENCIL_PROXY_ITERATOR_TRAITS_;
    typedef typename STENCIL_PROXY_ITERATOR_TRAITS_::ITERATOR_FACADE_ ITERATOR_FACADE_;
    typedef typename STENCIL_PROXY_ITERATOR_TRAITS_::INDEX_REFERENCE_TYPE INDEX_REFERENCE_TYPE;
    typedef typename STENCIL_PROXY_ITERATOR_TRAITS_::VALUE_REFERENCE_TYPE VALUE_REFERENCE_TYPE;
public:
    typedef typename ITERATOR_FACADE_::value_type value_type;
    typedef typename ITERATOR_FACADE_::reference reference;
    typedef typename ITERATOR_FACADE_::difference_type difference_type;

    typedef T_INDEX_ITERATOR INDEX_ITERATOR_TYPE;
    typedef T_VALUE_ITERATOR VALUE_ITERATOR_TYPE;

    typedef typename STENCIL_PROXY_ITERATOR_TRAITS_::INDEX_TYPE INDEX_TYPE;
    typedef typename STENCIL_PROXY_ITERATOR_TRAITS_::VALUE_TYPE VALUE_TYPE;

    STENCIL_PROXY_ITERATOR();
    STENCIL_PROXY_ITERATOR(const T_INDEX_ITERATOR& index_it, const T_VALUE_ITERATOR& value_it);

    template<class,class> friend class STENCIL_PROXY_ITERATOR;
    template< class T_INDEX_ITERATOR2, class T_VALUE_ITERATOR2 >
    STENCIL_PROXY_ITERATOR(const STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR2, T_VALUE_ITERATOR2 >& other,
        typename boost::enable_if< boost::mpl::and_<
            boost::is_convertible< T_INDEX_ITERATOR2, T_INDEX_ITERATOR >,
            boost::is_convertible< T_VALUE_ITERATOR2, T_VALUE_ITERATOR >
        > >::type* = 0);

    INDEX_REFERENCE_TYPE Index() const;
    VALUE_REFERENCE_TYPE Value() const;

    bool Valid() const;
    void Next();
    void Prev();

private:
    INDEX_ITERATOR_TYPE m_index_it;
    VALUE_ITERATOR_TYPE m_value_it;

    friend class boost::iterator_core_access;

    reference dereference() const;
    template< class T_OTHER >
    bool equal(const T_OTHER& other) const;
    void increment();
    void decrement();
    void advance(difference_type n);
    template< class T_OTHER >
    difference_type distance_to(const T_OTHER& other) const;
};

//#####################################################################
//#####################################################################

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
inline
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
STENCIL_PROXY_ITERATOR()
{ }

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
inline
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
STENCIL_PROXY_ITERATOR(const T_INDEX_ITERATOR& index_it, const T_VALUE_ITERATOR& value_it)
    : m_index_it(index_it),
      m_value_it(value_it)
{ }

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
template< class T_INDEX_ITERATOR2, class T_VALUE_ITERATOR2 >
inline
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
STENCIL_PROXY_ITERATOR(const STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR2, T_VALUE_ITERATOR2 >& other,
    typename boost::enable_if< boost::mpl::and_<
        boost::is_convertible< T_INDEX_ITERATOR2, T_INDEX_ITERATOR >,
        boost::is_convertible< T_VALUE_ITERATOR2, T_VALUE_ITERATOR >
    > >::type*)
    : m_index_it(other.m_index_it),
      m_value_it(other.m_value_it)
{ }

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
inline typename STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::INDEX_REFERENCE_TYPE
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
Index() const
{ return *m_index_it; }

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
inline typename STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::VALUE_REFERENCE_TYPE
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
Value() const
{ return *m_value_it; }

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
inline bool
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
Valid() const
{ return m_index_it.Valid(); }

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
inline void
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
Next()
{ increment(); }

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
inline void
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
Prev()
{ decrement(); }

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
inline typename STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::reference
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
dereference() const
{ return reference(*m_index_it, *m_value_it); }

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
template< class T_OTHER >
inline bool
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
equal(const T_OTHER& other) const
{
    assert((m_index_it == other.m_index_it) == (m_value_it == other.m_value_it));
    return m_value_it == other.m_value_it;
}

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
inline void
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
increment()
{
    ++m_index_it;
    ++m_value_it;
}

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
inline void
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
decrement()
{
    --m_index_it;
    --m_value_it;
}

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
inline void
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
advance(difference_type n)
{
    m_index_it += n;
    m_value_it += n;
}

template< class T_INDEX_ITERATOR, class T_VALUE_ITERATOR >
template< class T_OTHER >
inline typename STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::difference_type
STENCIL_PROXY_ITERATOR< T_INDEX_ITERATOR, T_VALUE_ITERATOR >::
distance_to(const T_OTHER& other) const
{
    assert(other.m_index_it - m_index_it == other.m_value_it - m_value_it);
    return static_cast< difference_type >(other.m_value_it - m_value_it);
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_STENCIL_PROXY_ITERATOR_HPP
