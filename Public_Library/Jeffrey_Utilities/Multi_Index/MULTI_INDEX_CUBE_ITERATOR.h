//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// struct MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_CUBE_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_CUBE_ITERATOR_HPP

#include <cassert>

#include <iterator>

#include <boost/mpl/assert.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/ITERATOR_FACADE.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_FWD.h>

#define PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX

namespace PhysBAM
{

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ = -MIN_OFFSET_ >
class MULTI_INDEX_CUBE_ITERATOR;

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
struct MULTI_INDEX_CUBE_ITERATOR_TRAITS
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef ITERATOR_FACADE<
        MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >, // Derived
        const MULTI_INDEX_TYPE,                                   // Value
        std::random_access_iterator_tag,                          // CategoryOrTraversal
        MULTI_INDEX_TYPE,                                         // Reference
        int                                                       // Difference
    > ITERATOR_FACADE_;
};

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
class MULTI_INDEX_CUBE_ITERATOR
    : public MULTI_INDEX_CUBE_ITERATOR_TRAITS< D, MIN_OFFSET_, MAX_OFFSET_ >::ITERATOR_FACADE_
{
    typedef MULTI_INDEX_CUBE_ITERATOR_TRAITS< D, MIN_OFFSET_, MAX_OFFSET_ > MULTI_INDEX_CUBE_ITERATOR_TRAITS_;
    typedef typename MULTI_INDEX_CUBE_ITERATOR_TRAITS_::ITERATOR_FACADE_ ITERATOR_FACADE_;
public:
    BOOST_MPL_ASSERT_RELATION( MIN_OFFSET_, <=, MAX_OFFSET_ );

    typedef typename ITERATOR_FACADE_::value_type value_type;
    typedef typename ITERATOR_FACADE_::reference reference;
    typedef typename ITERATOR_FACADE_::difference_type difference_type;

    typedef typename MULTI_INDEX_CUBE_ITERATOR_TRAITS_::MULTI_INDEX_TYPE MULTI_INDEX_TYPE;
    typedef MULTI_INDEX_CUBE< D, MIN_OFFSET_, MAX_OFFSET_ > MULTI_INDEX_CUBE_TYPE;

    static const int DIMENSION = D;
    static const int MIN_OFFSET = MIN_OFFSET_;
    static const int MAX_OFFSET = MAX_OFFSET_;
    static const int WIDTH = MULTI_INDEX_CUBE_TYPE::WIDTH;

    MULTI_INDEX_CUBE_ITERATOR();
    explicit MULTI_INDEX_CUBE_ITERATOR(const MULTI_INDEX_CUBE_TYPE& multi_index_cube);
    MULTI_INDEX_CUBE_ITERATOR(const MULTI_INDEX_CUBE_TYPE& multi_index_cube, BEGIN_TAG);
    MULTI_INDEX_CUBE_ITERATOR(const MULTI_INDEX_CUBE_TYPE& multi_index_cube, END_TAG);
    explicit MULTI_INDEX_CUBE_ITERATOR(const MULTI_INDEX_TYPE& base_multi_index);
    MULTI_INDEX_CUBE_ITERATOR(const MULTI_INDEX_TYPE& base_multi_index, BEGIN_TAG);
    MULTI_INDEX_CUBE_ITERATOR(const MULTI_INDEX_TYPE& base_multi_index, END_TAG);

    const MULTI_INDEX_CUBE_TYPE& Multi_Index_Cube() const;
    MULTI_INDEX_TYPE Multi_Index() const;
    MULTI_INDEX_TYPE Multi_Offset() const;
    int Linear_Index() const;

    bool Valid() const;
    void Next();
    void Prev();

private:
    MULTI_INDEX_CUBE_TYPE m_multi_index_cube;
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    int m_linear_index;
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    MULTI_INDEX_TYPE m_multi_index;
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX

    friend class boost::iterator_core_access;

    reference dereference() const;
    bool equal(const MULTI_INDEX_CUBE_ITERATOR& other) const;
    void increment();
    void decrement();
    void advance(difference_type n);
    difference_type distance_to(const MULTI_INDEX_CUBE_ITERATOR& other) const;
};

//##################################################################### 
//##################################################################### 

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
MULTI_INDEX_CUBE_ITERATOR()
{ }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
MULTI_INDEX_CUBE_ITERATOR(const MULTI_INDEX_CUBE_TYPE& multi_index_cube)
    : m_multi_index_cube(multi_index_cube),
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
      m_linear_index(1)
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
      m_multi_index(multi_index_cube.base_multi_index + MIN_OFFSET)
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{}

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
MULTI_INDEX_CUBE_ITERATOR(const MULTI_INDEX_CUBE_TYPE& multi_index_cube, BEGIN_TAG)
    : m_multi_index_cube(multi_index_cube),
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
      m_linear_index(1)
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
      m_multi_index(multi_index_cube.base_multi_index + MIN_OFFSET)
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{}

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
MULTI_INDEX_CUBE_ITERATOR(const MULTI_INDEX_CUBE_TYPE& multi_index_cube, END_TAG)
    : m_multi_index_cube(multi_index_cube),
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
      m_linear_index(1 + MULTI_INDEX_CUBE_TYPE::SIZE)
{ }
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
      m_multi_index(multi_index_cube.base_multi_index + MIN_OFFSET)
{ m_multi_index[1] = 1 + multi_index_cube.base_multi_index[1] + MAX_OFFSET; }
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
MULTI_INDEX_CUBE_ITERATOR(const MULTI_INDEX_TYPE& base_multi_index)
    : m_multi_index_cube(base_multi_index),
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
      m_linear_index(1)
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
      m_multi_index(base_multi_index + MIN_OFFSET)
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{}

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
MULTI_INDEX_CUBE_ITERATOR(const MULTI_INDEX_TYPE& base_multi_index, BEGIN_TAG)
    : m_multi_index_cube(base_multi_index),
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
      m_linear_index(1)
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
      m_multi_index(base_multi_index + MIN_OFFSET)
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{}

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
MULTI_INDEX_CUBE_ITERATOR(const MULTI_INDEX_TYPE& base_multi_index, END_TAG)
    : m_multi_index_cube(base_multi_index),
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
      m_linear_index(1 + MULTI_INDEX_CUBE_TYPE::SIZE)
{ }
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
      m_multi_index(base_multi_index + MIN_OFFSET)
{ m_multi_index[1] = 1 + base_multi_index[1] + MAX_OFFSET; }
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::MULTI_INDEX_CUBE_TYPE const &
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
Multi_Index_Cube() const
{ return m_multi_index_cube; }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::MULTI_INDEX_TYPE
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
Multi_Index() const
{ return dereference(); }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::MULTI_INDEX_TYPE
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
Multi_Offset() const
{ return m_multi_index_cube.Multi_Offset(dereference()); }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline int
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
Linear_Index() const
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{ return m_linear_index; }
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{ return m_multi_index_cube.Linear_Index(m_multi_index); }
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline bool
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
Valid() const
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{ return m_linear_index != 1 + MULTI_INDEX_CUBE_TYPE::SIZE; }
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{ return m_multi_index[1] != 1 + m_multi_index_cube.base_multi_index[1] + MAX_OFFSET; }
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline void
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
Next()
{ increment(); }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline void
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
Prev()
{ decrement(); }

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::reference
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
dereference() const
{
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    assert(1 <= m_linear_index && m_linear_index <= MULTI_INDEX_CUBE_TYPE::SIZE);
    return m_multi_index_cube.Multi_Index(m_linear_index);
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    assert(m_multi_index_cube.base_multi_index[1] + MIN_OFFSET <= m_multi_index[1] &&
           m_multi_index[1] <= m_multi_index_cube.base_multi_index[1] + MAX_OFFSET);
    return m_multi_index;
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
}

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline bool
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
equal(const MULTI_INDEX_CUBE_ITERATOR& other) const
{
    assert(m_multi_index_cube == other.m_multi_index_cube);
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    return m_linear_index == other.m_linear_index;
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    return m_multi_index == other.m_multi_index;
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
}

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline void
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
increment()
{
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    assert(m_linear_index <= MULTI_INDEX_CUBE_TYPE::SIZE);
    ++m_linear_index;
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    assert(m_multi_index[1] <= m_multi_index_cube.base_multi_index[1] + MAX_OFFSET);
    for(int d = D; d >= 2; --d) {
        if(++m_multi_index[d] <= m_multi_index_cube.base_multi_index[d] + MAX_OFFSET)
            return;
        m_multi_index[d] = m_multi_index_cube.base_multi_index[d] + MIN_OFFSET;
    }
    ++m_multi_index[1];
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
}

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline void
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
decrement()
{
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    --m_linear_index;
    assert(1 <= m_linear_index);
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    for(int d = D; d >= 2; --d) {
        if(--m_multi_index[d] >= m_multi_index_cube.base_multi_index[d] + MIN_OFFSET)
            return;
        m_multi_index[d] = m_multi_index_cube.base_multi_index[d] + MAX_OFFSET;
    }
    --m_multi_index[1];
    assert(m_multi_index_cube.base_multi_index[1] + MIN_OFFSET <= m_multi_index[1]);
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
}

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline void
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
advance(difference_type n)
{
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    m_linear_index += n;
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    if(n > 0) {
        for(int d = D; d >= 2; --d) {
            m_multi_index[d] += n % WIDTH;
            n /= WIDTH;
            if(m_multi_index[d] > m_multi_index_cube.base_multi_index[d] + MAX_OFFSET) {
                m_multi_index[d] -= WIDTH;
                ++n;
            }
        }
        m_multi_index[1] += n;
    }
    else if(n < 0) {
        n = -n;
        for(int d = D; d >= 2; --d) {
            m_multi_index[d] -= n % WIDTH;
            n /= WIDTH;
            if(m_multi_index[d] < m_multi_index_cube.base_multi_index[d] + MIN_OFFSET) {
                m_multi_index[d] += WIDTH;
                ++n;
            }
        }
        m_multi_index[1] -= n;
    }
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
}

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ >
inline typename MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::difference_type
MULTI_INDEX_CUBE_ITERATOR< D, MIN_OFFSET_, MAX_OFFSET_ >::
distance_to(const MULTI_INDEX_CUBE_ITERATOR& other) const
{
    assert(m_multi_index_cube == other.m_multi_index_cube);
#ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    return other.m_linear_index - m_linear_index;
#else // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    int n = other.m_multi_index[1] - m_multi_index[1];
    for(int d = 2; d <= D; ++d) {
        n *= WIDTH;
        n += other.m_multi_index[d] - m_multi_index[d];
    }
    return n;
#endif // #ifdef PHYSBAM_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_CUBE_ITERATOR_HPP
