//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// struct MULTI_INDEX_BOX_ITERATOR<D>
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_BOX_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_BOX_ITERATOR_HPP

#include <cassert>

#include <iterator>

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/ITERATOR_FACADE.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOX.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_FWD.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#define PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX

namespace PhysBAM
{

template< int D >
class MULTI_INDEX_BOX_ITERATOR;

template< int D >
struct MULTI_INDEX_BOX_ITERATOR_TRAITS
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef ITERATOR_FACADE<
        MULTI_INDEX_BOX_ITERATOR<D>,     // Derived
        const MULTI_INDEX_TYPE,          // Value
        std::random_access_iterator_tag, // CategoryOrTraversal
        MULTI_INDEX_TYPE,                // Reference
        int                              // Difference
    > ITERATOR_FACADE_;
};

template< int D >
class MULTI_INDEX_BOX_ITERATOR
    : public MULTI_INDEX_BOX_ITERATOR_TRAITS<D>::ITERATOR_FACADE_
{
    typedef MULTI_INDEX_BOX_ITERATOR_TRAITS<D> MULTI_INDEX_BOX_ITERATOR_TRAITS_;
    typedef typename MULTI_INDEX_BOX_ITERATOR_TRAITS_::ITERATOR_FACADE_ ITERATOR_FACADE_;
public:
    typedef typename ITERATOR_FACADE_::value_type value_type;
    typedef typename ITERATOR_FACADE_::reference reference;
    typedef typename ITERATOR_FACADE_::difference_type difference_type;

    typedef typename MULTI_INDEX_BOX_ITERATOR_TRAITS_::MULTI_INDEX_TYPE MULTI_INDEX_TYPE;
    typedef MULTI_INDEX_BOX<D> MULTI_INDEX_BOX_TYPE;

    static const int DIMENSION = D;

    MULTI_INDEX_BOX_ITERATOR();
    explicit MULTI_INDEX_BOX_ITERATOR(const MULTI_INDEX_BOX_TYPE& multi_index_box);
    MULTI_INDEX_BOX_ITERATOR(const MULTI_INDEX_BOX_TYPE& multi_index_box, BEGIN_TAG);
    MULTI_INDEX_BOX_ITERATOR(const MULTI_INDEX_BOX_TYPE& multi_index_box, END_TAG);

    const MULTI_INDEX_BOX_TYPE& Multi_Index_Box() const;
    MULTI_INDEX_TYPE Multi_Index() const;
    int Linear_Index() const;

    bool Valid() const;
    void Next();
    void Prev();

private:
    MULTI_INDEX_BOX_TYPE m_multi_index_box;
#ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    int m_linear_index;
#else // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    MULTI_INDEX_TYPE m_multi_index;
#endif // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX

    friend class boost::iterator_core_access;

    reference dereference() const;
    bool equal(const MULTI_INDEX_BOX_ITERATOR& other) const;
    void increment();
    void decrement();
    void advance(difference_type n);
    difference_type distance_to(const MULTI_INDEX_BOX_ITERATOR& other) const;
};

//##################################################################### 
//##################################################################### 

template< int D >
inline
MULTI_INDEX_BOX_ITERATOR<D>::
MULTI_INDEX_BOX_ITERATOR()
{ }

template< int D >
inline
MULTI_INDEX_BOX_ITERATOR<D>::
MULTI_INDEX_BOX_ITERATOR(const MULTI_INDEX_BOX_TYPE& multi_index_box)
    : m_multi_index_box(multi_index_box),
#ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
      m_linear_index(1)
#else // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
      m_multi_index(multi_index_box.min_multi_index)
#endif // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
{}

template< int D >
inline
MULTI_INDEX_BOX_ITERATOR<D>::
MULTI_INDEX_BOX_ITERATOR(const MULTI_INDEX_BOX_TYPE& multi_index_box, BEGIN_TAG)
    : m_multi_index_box(multi_index_box),
#ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
      m_linear_index(1)
#else // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
      m_multi_index(multi_index_box.min_multi_index)
#endif // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
{}

template< int D >
inline
MULTI_INDEX_BOX_ITERATOR<D>::
MULTI_INDEX_BOX_ITERATOR(const MULTI_INDEX_BOX_TYPE& multi_index_box, END_TAG)
    : m_multi_index_box(multi_index_box),
#ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
      m_linear_index(1 + multi_index_box.Size())
{ }
#else // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
      m_multi_index(multi_index_box.min_multi_index)
{ m_multi_index[1] = 1 + multi_index_box.max_multi_index[1]; }
#endif // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX

template< int D >
inline typename MULTI_INDEX_BOX_ITERATOR<D>::MULTI_INDEX_BOX_TYPE const &
MULTI_INDEX_BOX_ITERATOR<D>::
Multi_Index_Box() const
{ return m_multi_index_box; }

template< int D >
inline typename MULTI_INDEX_BOX_ITERATOR<D>::MULTI_INDEX_TYPE
MULTI_INDEX_BOX_ITERATOR<D>::
Multi_Index() const
{ return dereference(); }

template< int D >
inline int
MULTI_INDEX_BOX_ITERATOR<D>::
Linear_Index() const
#ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
{ return m_linear_index; }
#else // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
{ return m_multi_index_box.Linear_Index(m_multi_index); }
#endif // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX

template< int D >
inline bool
MULTI_INDEX_BOX_ITERATOR<D>::
Valid() const
#ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
{ return m_linear_index != 1 + m_multi_index_box.Size(); }
#else // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
{ return m_multi_index[1] != 1 + m_multi_index_box.max_multi_index[1]; }
#endif // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX

template< int D >
inline void
MULTI_INDEX_BOX_ITERATOR<D>::
Next()
{ increment(); }

template< int D >
inline void
MULTI_INDEX_BOX_ITERATOR<D>::
Prev()
{ decrement(); }

template< int D >
inline typename MULTI_INDEX_BOX_ITERATOR<D>::reference
MULTI_INDEX_BOX_ITERATOR<D>::
dereference() const
{
#ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    assert(1 <= m_linear_index && m_linear_index <= m_multi_index_box.Size());
    return m_multi_index_box.Multi_Index(m_linear_index);
#else // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    assert(m_multi_index_box.min_multi_index[1] <= m_multi_index[1] &&
           m_multi_index[1] <= m_multi_index_box.max_multi_index[1]);
    return m_multi_index;
#endif // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
}

template< int D >
inline bool
MULTI_INDEX_BOX_ITERATOR<D>::
equal(const MULTI_INDEX_BOX_ITERATOR& other) const
{
    assert(m_multi_index_box == other.m_multi_index_box);
#ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    return m_linear_index == other.m_linear_index;
#else // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    return m_multi_index == other.m_multi_index;
#endif // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
}

template< int D >
inline void
MULTI_INDEX_BOX_ITERATOR<D>::
increment()
{
#ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    assert(m_linear_index <= m_multi_index_box.Size());
    ++m_linear_index;
#else // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    assert(m_multi_index[1] <= m_multi_index_box.max_multi_index[1]);
    for(int d = D; d >= 2; --d) {
        if(++m_multi_index[d] <= m_multi_index_box.max_multi_index[d])
            return;
        m_multi_index[d] = m_multi_index_box.min_multi_index[d];
    }
    ++m_multi_index[1];
#endif // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
}

template< int D >
inline void
MULTI_INDEX_BOX_ITERATOR<D>::
decrement()
{
#ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    --m_linear_index;
    assert(1 <= m_linear_index);
#else // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    for(int d = D; d >= 2; --d) {
        if(--m_multi_index[d] >= m_multi_index_box.min_multi_index[d])
            return;
        m_multi_index[d] = m_multi_index_box.max_multi_index[d];
    }
    --m_multi_index[1];
    assert(m_multi_index_box.min_multi_index[1] <= m_multi_index[1]);
#endif // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
}

template< int D >
inline void
MULTI_INDEX_BOX_ITERATOR<D>::
advance(difference_type n)
{
#ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    m_linear_index += n;
#else // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    if(n > 0) {
        for(int d = D; d >= 2; --d) {
            const int w = m_multi_index_box.widths[d];
            m_multi_index[d] += n % w;
            n /= w;
            if(m_multi_index[d] > m_multi_index_box.max_multi_index[d]) {
                m_multi_index[d] -= w;
                ++n;
            }
        }
        m_multi_index[1] += n;
    }
    else if(n < 0) {
        n = -n;
        for(int d = D; d >= 2; --d) {
            const int w = m_multi_index_box.widths[d];
            m_multi_index[d] -= n % w;
            n /= w;
            if(m_multi_index[d] < m_multi_index_box.min_multi_index[d]) {
                m_multi_index[d] += w;
                ++n;
            }
        }
        m_multi_index[1] -= n;
    }
#endif // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
}

template< int D >
inline typename MULTI_INDEX_BOX_ITERATOR<D>::difference_type
MULTI_INDEX_BOX_ITERATOR<D>::
distance_to(const MULTI_INDEX_BOX_ITERATOR& other) const
{
    assert(m_multi_index_box == other.m_multi_index_box);
#ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    return other.m_linear_index - m_linear_index;
#else // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
    int n = other.m_multi_index[1] - m_multi_index[1];
    for(int d = 2; d <= D; ++d) {
        n *= m_multi_index_box.widths[d];
        n += other.m_multi_index[d] - m_multi_index[d];
    }
    return n;
#endif // #ifdef PHYSBAM_MULTI_INDEX_BOX_ITERATOR_USE_LINEAR_INDEX
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_BOX_ITERATOR_HPP
