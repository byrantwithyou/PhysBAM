//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_IDENTITY_STENCIL_PROXY_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_IDENTITY_STENCIL_PROXY_ITERATOR_HPP

#include <cassert>

#include <iterator>

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/ITERATOR_FACADE.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>

namespace PhysBAM
{

template< class T_INDEX, class T >
class IDENTITY_STENCIL_PROXY_ITERATOR;

template< class T_INDEX, class T >
struct IDENTITY_STENCIL_PROXY_ITERATOR_TRAITS
{
    typedef INDEX_VALUE< T_INDEX, T > INDEX_VALUE_TYPE;
    typedef ITERATOR_FACADE<
        IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >, // Derived
        const INDEX_VALUE_TYPE,                        // Value
        std::random_access_iterator_tag,               // CategoryOrTraversal
        INDEX_VALUE_TYPE,                              // Reference
        int                                            // Difference
    > ITERATOR_FACADE_;
};

template< class T_INDEX, class T >
class IDENTITY_STENCIL_PROXY_ITERATOR
    : public IDENTITY_STENCIL_PROXY_ITERATOR_TRAITS< T_INDEX, T >::ITERATOR_FACADE_
{
    typedef IDENTITY_STENCIL_PROXY_ITERATOR_TRAITS< T_INDEX, T > IDENTITY_STENCIL_PROXY_ITERATOR_TRAITS_;
    typedef typename IDENTITY_STENCIL_PROXY_ITERATOR_TRAITS_::ITERATOR_FACADE_ ITERATOR_FACADE_;
public:
    typedef typename ITERATOR_FACADE_::value_type value_type;
    typedef typename ITERATOR_FACADE_::reference reference;
    typedef typename ITERATOR_FACADE_::difference_type difference_type;

    typedef T_INDEX INDEX_TYPE;
    typedef T SCALAR_TYPE;

    IDENTITY_STENCIL_PROXY_ITERATOR();
    explicit IDENTITY_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& index);
    IDENTITY_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& index, BEGIN_TAG);
    IDENTITY_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& index, END_TAG);

    const T_INDEX& Index() const;

    bool Valid() const;
    void Next();
    void Prev();

private:
    INDEX_TYPE m_index;
    bool m_valid;

    friend class boost::iterator_core_access;

    reference dereference() const;
    bool equal(const IDENTITY_STENCIL_PROXY_ITERATOR& other) const;
    void increment();
    void decrement();
    void advance(difference_type n);
    difference_type distance_to(const IDENTITY_STENCIL_PROXY_ITERATOR& other) const;
};

//##################################################################### 
//##################################################################### 

template< class T_INDEX, class T >
inline
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
IDENTITY_STENCIL_PROXY_ITERATOR()
{ }

template< class T_INDEX, class T >
inline
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
IDENTITY_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& index)
    : m_index(index),
      m_valid(true)
{ }

template< class T_INDEX, class T >
inline
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
IDENTITY_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& index, BEGIN_TAG)
    : m_index(index),
      m_valid(true)
{}

template< class T_INDEX, class T >
inline
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
IDENTITY_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& index, END_TAG)
    : m_index(index),
      m_valid(false)
{}

template< class T_INDEX, class T >
inline const T_INDEX&
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
Index() const
{ return m_index; }

template< class T_INDEX, class T >
inline bool
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
Valid() const
{ return m_valid; }

template< class T_INDEX, class T >
inline void
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
Next()
{ increment(); }

template< class T_INDEX, class T >
inline void
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
Prev()
{ decrement(); }

template< class T_INDEX, class T >
inline typename IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::reference
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
dereference() const
{
    assert(m_valid);
    return reference(m_index, static_cast<T>(1));
}

template< class T_INDEX, class T >
inline bool
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
equal(const IDENTITY_STENCIL_PROXY_ITERATOR& other) const
{
    assert(m_index == other.m_index);
    return m_valid == other.m_valid;
}

template< class T_INDEX, class T >
inline void
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
increment()
{
    assert(m_valid);
    m_valid = false;
}

template< class T_INDEX, class T >
inline void
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
decrement()
{
    assert(!m_valid);
    m_valid = true;
}

template< class T_INDEX, class T >
inline void
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
advance(difference_type n)
{
    if(n != 0) {
        assert(m_valid && n == 1 || !m_valid && n == -1);
        m_valid = !m_valid;
    }
}

template< class T_INDEX, class T >
inline typename IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::difference_type
IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T >::
distance_to(const IDENTITY_STENCIL_PROXY_ITERATOR& other) const
{
    assert(m_index == other.m_index);
    return static_cast< int >(m_valid) - static_cast< int >(other.m_valid);
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_IDENTITY_STENCIL_PROXY_ITERATOR_HPP
