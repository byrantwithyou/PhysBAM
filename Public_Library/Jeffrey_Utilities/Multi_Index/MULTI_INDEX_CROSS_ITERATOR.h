//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// struct MULTI_INDEX_CROSS_ITERATOR<D>
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_CROSS_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_CROSS_ITERATOR_HPP

#include <cassert>
#include <cstdlib>

#include <iterator>

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/ITERATOR_FACADE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< int D >
class MULTI_INDEX_CROSS_ITERATOR;

template< int D >
struct MULTI_INDEX_CROSS_ITERATOR_TRAITS
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef ITERATOR_FACADE<
        MULTI_INDEX_CROSS_ITERATOR<D>,   // Derived
        const MULTI_INDEX_TYPE,          // Value
        std::random_access_iterator_tag, // CategoryOrTraversal
        MULTI_INDEX_TYPE,                // Reference
        int                              // Difference
    > ITERATOR_FACADE_;
};

template< int D >
class MULTI_INDEX_CROSS_ITERATOR
    : public MULTI_INDEX_CROSS_ITERATOR_TRAITS<D>::ITERATOR_FACADE_
{
    typedef MULTI_INDEX_CROSS_ITERATOR_TRAITS<D> MULTI_INDEX_CROSS_ITERATOR_TRAITS_;
    typedef typename MULTI_INDEX_CROSS_ITERATOR_TRAITS_::ITERATOR_FACADE_ ITERATOR_FACADE_;
public:
    typedef typename ITERATOR_FACADE_::value_type value_type;
    typedef typename ITERATOR_FACADE_::reference reference;
    typedef typename ITERATOR_FACADE_::difference_type difference_type;

    typedef typename MULTI_INDEX_CROSS_ITERATOR_TRAITS_::MULTI_INDEX_TYPE MULTI_INDEX_TYPE;

    MULTI_INDEX_CROSS_ITERATOR();
    explicit MULTI_INDEX_CROSS_ITERATOR(const MULTI_INDEX_TYPE& base_multi_index);
    MULTI_INDEX_CROSS_ITERATOR(const MULTI_INDEX_TYPE& base_multi_index, BEGIN_TAG);
    MULTI_INDEX_CROSS_ITERATOR(const MULTI_INDEX_TYPE& base_multi_index, END_TAG);

    MULTI_INDEX_TYPE Multi_Index() const;
    int Linear_Index() const;

    bool Valid() const;
    void Next();
    void Prev();

private:
    MULTI_INDEX_TYPE m_base_multi_index;
    int m_index;

    friend class boost::iterator_core_access;

    reference dereference() const;
    bool equal(const MULTI_INDEX_CROSS_ITERATOR& other) const;
    void increment();
    void decrement();
    void advance(const difference_type n);
    difference_type distance_to(const MULTI_INDEX_CROSS_ITERATOR& other) const;
};

//##################################################################### 
//##################################################################### 

template< int D >
inline
MULTI_INDEX_CROSS_ITERATOR<D>::
MULTI_INDEX_CROSS_ITERATOR()
{ }

template< int D >
inline
MULTI_INDEX_CROSS_ITERATOR<D>::
MULTI_INDEX_CROSS_ITERATOR(const MULTI_INDEX_TYPE& base_multi_index)
    : m_base_multi_index(base_multi_index),
      m_index(0)
{ }

template< int D >
inline
MULTI_INDEX_CROSS_ITERATOR<D>::
MULTI_INDEX_CROSS_ITERATOR(const MULTI_INDEX_TYPE& base_multi_index, BEGIN_TAG)
    : m_base_multi_index(base_multi_index),
      m_index(0)
{ }

template< int D >
inline
MULTI_INDEX_CROSS_ITERATOR<D>::
MULTI_INDEX_CROSS_ITERATOR(const MULTI_INDEX_TYPE& base_multi_index, END_TAG)
    : m_base_multi_index(base_multi_index),
      m_index(2*D + 1)
{ }

template< int D >
inline typename MULTI_INDEX_CROSS_ITERATOR<D>::MULTI_INDEX_TYPE
MULTI_INDEX_CROSS_ITERATOR<D>::
Multi_Index() const
{ return dereference(); }

template< int D >
inline int
MULTI_INDEX_CROSS_ITERATOR<D>::
Linear_Index() const
{ return m_index + 1; }

template< int D >
inline bool
MULTI_INDEX_CROSS_ITERATOR<D>::
Valid() const
{ return m_index != 2*D + 1; }

template< int D >
inline void
MULTI_INDEX_CROSS_ITERATOR<D>::
Next()
{ increment(); }

template< int D >
inline void
MULTI_INDEX_CROSS_ITERATOR<D>::
Prev()
{ decrement(); }

template< int D >
inline typename MULTI_INDEX_CROSS_ITERATOR<D>::reference
MULTI_INDEX_CROSS_ITERATOR<D>::
dereference() const
{
    assert(0 <= m_index);
    assert(m_index <= 2*D);
    MULTI_INDEX_TYPE multi_index = m_base_multi_index;
    if(const int j = std::abs(m_index - D))
        multi_index[D - j + 1] += 2*(D < m_index) - 1;
    return multi_index;
}

template< int D >
inline bool
MULTI_INDEX_CROSS_ITERATOR<D>::
equal(const MULTI_INDEX_CROSS_ITERATOR& other) const
{
    assert(m_base_multi_index == other.m_base_multi_index);
    return m_index == other.m_index;
}

template< int D >
inline void
MULTI_INDEX_CROSS_ITERATOR<D>::
increment()
{
    assert(m_index <= 2*D);
    ++m_index;
}

template< int D >
inline void
MULTI_INDEX_CROSS_ITERATOR<D>::
decrement()
{
    --m_index;
    assert(0 <= m_index);
}

template< int D >
inline void
MULTI_INDEX_CROSS_ITERATOR<D>::
advance(const difference_type n)
{ m_index += n; }

template< int D >
inline typename MULTI_INDEX_CROSS_ITERATOR<D>::difference_type
MULTI_INDEX_CROSS_ITERATOR<D>::
distance_to(const MULTI_INDEX_CROSS_ITERATOR& other) const
{
    assert(m_base_multi_index == other.m_base_multi_index);
    return other.m_index - m_index;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_CROSS_ITERATOR_HPP
