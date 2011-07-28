//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_STATIC_MULTI_INDEX_CUBE_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_STATIC_MULTI_INDEX_CUBE_ITERATOR_HPP

#include <cassert>

#include <iterator>

#include <boost/mpl/assert.hpp>

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/ITERATOR_FACADE.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_FWD.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#define PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX

namespace PhysBAM
{

template< int D, int MIN_INDEX_, int MAX_INDEX_ = -MIN_INDEX_ >
class STATIC_MULTI_INDEX_CUBE_ITERATOR;

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
struct STATIC_MULTI_INDEX_CUBE_ITERATOR_TRAITS
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;
    typedef ITERATOR_FACADE<
        STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >, // Derived
        const MULTI_INDEX_TYPE,                                        // Value
        std::random_access_iterator_tag,                               // CategoryOrTraversal
        MULTI_INDEX_TYPE,                                              // Reference
        int                                                            // Difference
    > ITERATOR_FACADE_;
};

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
class STATIC_MULTI_INDEX_CUBE_ITERATOR
    : public STATIC_MULTI_INDEX_CUBE_ITERATOR_TRAITS< D, MIN_INDEX_, MAX_INDEX_ >::ITERATOR_FACADE_
{
    typedef STATIC_MULTI_INDEX_CUBE_ITERATOR_TRAITS< D, MIN_INDEX_, MAX_INDEX_ > STATIC_MULTI_INDEX_CUBE_ITERATOR_TRAITS_;
    typedef typename STATIC_MULTI_INDEX_CUBE_ITERATOR_TRAITS_::ITERATOR_FACADE_ ITERATOR_FACADE_;
public:
    BOOST_MPL_ASSERT_RELATION( MIN_INDEX_, <=, MAX_INDEX_ );

    typedef typename ITERATOR_FACADE_::value_type value_type;
    typedef typename ITERATOR_FACADE_::reference reference;
    typedef typename ITERATOR_FACADE_::difference_type difference_type;

    typedef typename STATIC_MULTI_INDEX_CUBE_ITERATOR_TRAITS_::MULTI_INDEX_TYPE MULTI_INDEX_TYPE;
    typedef STATIC_MULTI_INDEX_CUBE< D, MIN_INDEX_, MAX_INDEX_ > STATIC_MULTI_INDEX_CUBE_TYPE;

    static const int DIMENSION = D;
    static const int MIN_INDEX = MIN_INDEX_;
    static const int MAX_INDEX = MAX_INDEX_;
    static const int WIDTH = STATIC_MULTI_INDEX_CUBE_TYPE::WIDTH;

    STATIC_MULTI_INDEX_CUBE_ITERATOR();
    explicit STATIC_MULTI_INDEX_CUBE_ITERATOR(BEGIN_TAG);
    explicit STATIC_MULTI_INDEX_CUBE_ITERATOR(END_TAG);

    MULTI_INDEX_TYPE Multi_Index() const;
    int Linear_Index() const;

    bool Valid() const;
    void Next();
    void Prev();

private:
#ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    int m_linear_index;
#else // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    MULTI_INDEX_TYPE m_multi_index;
#endif // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX

    friend class boost::iterator_core_access;

    reference dereference() const;
    bool equal(const STATIC_MULTI_INDEX_CUBE_ITERATOR& other) const;
    void increment();
    void decrement();
    void advance(difference_type n);
    difference_type distance_to(const STATIC_MULTI_INDEX_CUBE_ITERATOR& other) const;
};

//##################################################################### 
//##################################################################### 

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
STATIC_MULTI_INDEX_CUBE_ITERATOR()
#ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    : m_linear_index(1)
#else // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    : m_multi_index(MIN_INDEX * MULTI_INDEX_TYPE::All_Ones_Vector())
#endif // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{ }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
STATIC_MULTI_INDEX_CUBE_ITERATOR(BEGIN_TAG)
#ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    : m_linear_index(1)
#else // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    : m_multi_index(MIN_INDEX * MULTI_INDEX_TYPE::All_Ones_Vector())
#endif // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{}

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
STATIC_MULTI_INDEX_CUBE_ITERATOR(END_TAG)
#ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    : m_linear_index(1 + STATIC_MULTI_INDEX_CUBE_TYPE::SIZE)
{ }
#else // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    : m_multi_index(MIN_INDEX * MULTI_INDEX_TYPE::All_Ones_Vector())
{ m_multi_index[1] = 1 + MAX_INDEX; }
#endif // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline typename STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::MULTI_INDEX_TYPE
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
Multi_Index() const
{ return dereference(); }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline int
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
Linear_Index() const
#ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{ return m_linear_index; }
#else // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{ return STATIC_MULTI_INDEX_CUBE_TYPE::Linear_Index(m_multi_index); }
#endif // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline bool
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
Valid() const
#ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{ return m_linear_index != 1 + STATIC_MULTI_INDEX_CUBE_TYPE::SIZE; }
#else // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
{ return m_multi_index[1] != 1 + MAX_INDEX; }
#endif // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline void
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
Next()
{ increment(); }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline void
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
Prev()
{ decrement(); }

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline typename STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::reference
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
dereference() const
{
#ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    assert(1 <= m_linear_index && m_linear_index <= STATIC_MULTI_INDEX_CUBE_TYPE::SIZE);
    return STATIC_MULTI_INDEX_CUBE_TYPE::Multi_Index(m_linear_index);
#else // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    assert(MIN_INDEX <= m_multi_index[1] && m_multi_index[1] <= MAX_INDEX);
    return m_multi_index;
#endif // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
}

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline bool
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
equal(const STATIC_MULTI_INDEX_CUBE_ITERATOR& other) const
{
#ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    return m_linear_index == other.m_linear_index;
#else // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    return m_multi_index == other.m_multi_index;
#endif // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
}

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline void
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
increment()
{
#ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    assert(m_linear_index <= STATIC_MULTI_INDEX_CUBE_TYPE::SIZE);
    ++m_linear_index;
#else // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    assert(m_multi_index[1] <= MAX_INDEX);
    for(int d = D; d >= 2; --d) {
        if(++m_multi_index[d] <= MAX_INDEX)
            return;
        m_multi_index[d] = MIN_INDEX;
    }
    ++m_multi_index[1];
#endif // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
}

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline void
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
decrement()
{
#ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    --m_linear_index;
    assert(1 <= m_linear_index);
#else // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    for(int d = D; d >= 2; --d) {
        if(--m_multi_index[d] >= MIN_INDEX)
            return;
        m_multi_index[d] = MAX_INDEX;
    }
    --m_multi_index[1];
    assert(MIN_INDEX <= m_multi_index[1]);
#endif // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
}

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline void
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
advance(difference_type n)
{
#ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    m_linear_index += n;
#else // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    if(n > 0) {
        for(int d = D; d >= 2; --d) {
            m_multi_index[d] += n % WIDTH;
            n /= WIDTH;
            if(m_multi_index[d] > MAX_INDEX) {
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
            if(m_multi_index[d] < MIN_INDEX) {
                m_multi_index[d] += WIDTH;
                ++n;
            }
        }
        m_multi_index[1] -= n;
    }
#endif // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
}

template< int D, int MIN_INDEX_, int MAX_INDEX_ >
inline typename STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::difference_type
STATIC_MULTI_INDEX_CUBE_ITERATOR< D, MIN_INDEX_, MAX_INDEX_ >::
distance_to(const STATIC_MULTI_INDEX_CUBE_ITERATOR& other) const
{
#ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    return other.m_linear_index - m_linear_index;
#else // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
    int n = other.m_multi_index[1] - m_multi_index[1];
    for(int d = 2; d <= D; ++d) {
        n *= WIDTH;
        n += other.m_multi_index[d] - m_multi_index[d];
    }
    return n;
#endif // #ifdef PHYSBAM_STATIC_MULTI_INDEX_CUBE_ITERATOR_USE_LINEAR_INDEX
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_STATIC_MULTI_INDEX_CUBE_ITERATOR_HPP
