//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_CUBE2_SIMPLEX_PARTITION_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_CUBE2_SIMPLEX_PARTITION_ITERATOR_HPP

#include <cassert>

#include <algorithm>
#include <iterator>

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/ITERATOR_FACADE.h>
#include <Jeffrey_Utilities/Math/STATIC_POW.h>
#include <Jeffrey_Utilities/Multi_Index/STATIC_MULTI_INDEX_CUBE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< int D >
class CUBE2_SIMPLEX_PARTITION_ITERATOR;

template< int D >
struct CUBE2_SIMPLEX_PARTITION_ITERATOR_TRAITS
{
    typedef ITERATOR_FACADE<
        CUBE2_SIMPLEX_PARTITION_ITERATOR<D>, // Derived
        const VECTOR< int, D+1 >,            // Value
        //std::random_access_iterator_tag,     // CategoryOrTraversal
        std::forward_iterator_tag,           // CategoryOrTraversal
        VECTOR< int, D+1 >,                  // Reference
        int                                  // Difference
    > ITERATOR_FACADE_;
};

template<>
class CUBE2_SIMPLEX_PARTITION_ITERATOR<1>
    : public CUBE2_SIMPLEX_PARTITION_ITERATOR_TRAITS<1>::ITERATOR_FACADE_
{
    typedef CUBE2_SIMPLEX_PARTITION_ITERATOR_TRAITS<1> CUBE2_SIMPLEX_PARTITION_ITERATOR_TRAITS_;
    typedef CUBE2_SIMPLEX_PARTITION_ITERATOR_TRAITS_::ITERATOR_FACADE_ ITERATOR_FACADE_;
public:
    typedef ITERATOR_FACADE_::value_type value_type;
    typedef ITERATOR_FACADE_::reference reference;
    typedef ITERATOR_FACADE_::difference_type difference_type;

    static const int DIMENSION = 1;
    static const int CENTER_INDEX = 2;
    static const int N_FACE = 2;
    static const int N_SIMPLEX = 2;

    CUBE2_SIMPLEX_PARTITION_ITERATOR();
    explicit CUBE2_SIMPLEX_PARTITION_ITERATOR(BEGIN_TAG);
    explicit CUBE2_SIMPLEX_PARTITION_ITERATOR(END_TAG);

    bool At_Begin() const;
    bool At_End() const;
    void Set_Begin();
    void Set_End();

    bool Valid() const;
    void Next();
    void Prev();

private:
    int m_face_index;

    friend class boost::iterator_core_access;

    reference dereference() const;
    bool equal(const CUBE2_SIMPLEX_PARTITION_ITERATOR& other) const;
    void increment();
    void decrement();
    //void advance(const difference_type n);
    //difference_type distance_to(const CUBE2_SIMPLEX_PARTITION_ITERATOR& other) const;
};

template< int D >
class CUBE2_SIMPLEX_PARTITION_ITERATOR
    : public CUBE2_SIMPLEX_PARTITION_ITERATOR_TRAITS<D>::ITERATOR_FACADE_
{
    typedef CUBE2_SIMPLEX_PARTITION_ITERATOR_TRAITS<D> CUBE2_SIMPLEX_PARTITION_ITERATOR_TRAITS_;
    typedef typename CUBE2_SIMPLEX_PARTITION_ITERATOR_TRAITS_::ITERATOR_FACADE_ ITERATOR_FACADE_;
public:
    typedef typename ITERATOR_FACADE_::value_type value_type;
    typedef typename ITERATOR_FACADE_::reference reference;
    typedef typename ITERATOR_FACADE_::difference_type difference_type;

    static const int DIMENSION = D;
    static const int CENTER_INDEX = (STATIC_POW_C<3,D>::value + 1) / 2;
    static const int N_FACE = 2 * D;
    static const int N_SIMPLEX = N_FACE * CUBE2_SIMPLEX_PARTITION_ITERATOR< D-1 >::N_SIMPLEX;

    CUBE2_SIMPLEX_PARTITION_ITERATOR();
    explicit CUBE2_SIMPLEX_PARTITION_ITERATOR(BEGIN_TAG);
    explicit CUBE2_SIMPLEX_PARTITION_ITERATOR(END_TAG);

    bool At_Begin() const;
    bool At_End() const;
    void Set_Begin();
    void Set_End();

    bool Valid() const;
    void Next();
    void Prev();

private:
    int m_face_index;
    CUBE2_SIMPLEX_PARTITION_ITERATOR< D-1 > m_sub_simplex_it;

    friend class boost::iterator_core_access;

    reference dereference() const;
    bool equal(const CUBE2_SIMPLEX_PARTITION_ITERATOR& other) const;
    void increment();
    void decrement();
    //void advance(const difference_type n);
    //difference_type distance_to(const CUBE2_SIMPLEX_PARTITION_ITERATOR& other) const;
};

//#####################################################################
//#####################################################################

inline
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
CUBE2_SIMPLEX_PARTITION_ITERATOR()
    : m_face_index(1)
{ }

inline
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
CUBE2_SIMPLEX_PARTITION_ITERATOR(BEGIN_TAG)
    : m_face_index(1)
{ }

inline
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
CUBE2_SIMPLEX_PARTITION_ITERATOR(END_TAG)
    : m_face_index(1 + N_FACE)
{ }

inline bool
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
At_Begin() const
{ return m_face_index == 1; }

inline bool
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
At_End() const
{ return m_face_index == 1 + N_FACE; }

inline void
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
Set_Begin()
{ m_face_index = 1; }

inline void
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
Set_End()
{ m_face_index = 1 + N_FACE; }

inline bool
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
Valid() const
{ return m_face_index <= N_FACE; }

inline void
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
Next()
{ increment(); }

inline void
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
Prev()
{ decrement(); }

inline CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::reference
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
dereference() const
{ return reference(m_face_index + 1, m_face_index); }

inline bool
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
equal(const CUBE2_SIMPLEX_PARTITION_ITERATOR& other) const
{ return m_face_index == other.m_face_index; }

inline void
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
increment()
{
    assert(m_face_index <= N_FACE);
    ++m_face_index;
}

inline void
CUBE2_SIMPLEX_PARTITION_ITERATOR<1>::
decrement()
{
    --m_face_index;
    assert(1 <= m_face_index);
}

//#####################################################################
//#####################################################################

template< int D >
inline
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
CUBE2_SIMPLEX_PARTITION_ITERATOR()
    : m_face_index(1)
{ }

template< int D >
inline
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
CUBE2_SIMPLEX_PARTITION_ITERATOR(BEGIN_TAG)
    : m_face_index(1),
      m_sub_simplex_it(BEGIN_TAG())
{ }

template< int D >
inline
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
CUBE2_SIMPLEX_PARTITION_ITERATOR(END_TAG)
    : m_face_index(1 + N_FACE),
      m_sub_simplex_it(BEGIN_TAG())
{ }

template< int D >
inline bool
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
At_Begin() const
{ return m_face_index == 1 && m_sub_simplex_it.At_Begin(); }

template< int D >
inline bool
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
At_End() const
{ return m_face_index == 1 + N_FACE; }

template< int D >
inline void
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
Set_Begin()
{
    m_face_index = 1;
    m_sub_simplex_it.Set_Begin();
}

template< int D >
inline void
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
Set_End()
{
    m_face_index = 1 + N_FACE;
    m_sub_simplex_it.Set_Begin();
}

template< int D >
inline bool
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
Valid() const
{ return m_face_index <= N_FACE; }

template< int D >
inline void
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
Next()
{ increment(); }

template< int D >
inline void
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
Prev()
{ decrement(); }

template< int D >
inline typename CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::reference
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
dereference() const
{
    const int axis_index = (m_face_index - 1) / 2;
    const int sign = 2 * ((m_face_index - 1) & 1) - 1;
    VECTOR<int,D> sub_simplex = *m_sub_simplex_it;
    for(int d = 1; d <= D; ++d) {
        const VECTOR< int, D-1 > sub_multi_index = STATIC_MULTI_INDEX_CUBE<D-1,-1,+1>::Multi_Index(sub_simplex[d]);
        const VECTOR<int,D> multi_index = sub_multi_index.Insert(sign, 1 + axis_index);
        sub_simplex[d] = STATIC_MULTI_INDEX_CUBE<D,-1,+1>::Linear_Index(multi_index);
    }
    // GCC workaround...???
    static const int CENTER_INDEX_ = CENTER_INDEX;
    assert(!sub_simplex.Contains(CENTER_INDEX_));
    VECTOR< int, D+1 > simplex = sub_simplex.Insert(CENTER_INDEX_, 1);
    if(((axis_index ^ (sign+1)/2) & 1) == 1)
        std::swap(simplex[D], simplex[D+1]);
    return simplex;
}

template< int D >
inline bool
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
equal(const CUBE2_SIMPLEX_PARTITION_ITERATOR& other) const
{ return m_face_index == other.m_face_index && m_sub_simplex_it == other.m_sub_simplex_it; }

template< int D >
inline void
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
increment()
{
    if((++m_sub_simplex_it).At_End()) {
        assert(m_face_index <= N_FACE);
        ++m_face_index;
        m_sub_simplex_it.Set_Begin();
    }
}

template< int D >
inline void
CUBE2_SIMPLEX_PARTITION_ITERATOR<D>::
decrement()
{
    if(m_sub_simplex_it.At_Begin()) {
        --m_face_index;
        assert(1 <= m_face_index);
        m_sub_simplex_it.Set_End();
    }
    --m_sub_simplex_it;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_CUBE2_SIMPLEX_PARTITION_ITERATOR_HPP
