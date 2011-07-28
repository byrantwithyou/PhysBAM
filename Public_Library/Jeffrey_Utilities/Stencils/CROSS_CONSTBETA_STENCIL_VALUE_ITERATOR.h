//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR_HPP

#include <cassert>
#include <cstdlib>

#include <iterator>

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/ITERATOR_FACADE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int D >
class CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR;

template< class T, int D >
struct CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR_TRAITS
{
    typedef ITERATOR_FACADE<
        CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>, // Derived
        const T,                                     // Value
        std::random_access_iterator_tag,             // CategoryOrTraversal
        T,                                           // Reference
        int                                          // Difference
    > ITERATOR_FACADE_;
};

template< class T, int D >
class CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR
    : public CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR_TRAITS<T,D>::ITERATOR_FACADE_
{
    typedef CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR_TRAITS<T,D> CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR_TRAITS_;
    typedef typename CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR_TRAITS_::ITERATOR_FACADE_ ITERATOR_FACADE_;
public:
    typedef typename ITERATOR_FACADE_::value_type value_type;
    typedef typename ITERATOR_FACADE_::reference reference;
    typedef typename ITERATOR_FACADE_::difference_type difference_type;

    typedef T SCALAR_TYPE;
    static const int DIMENSION = D;

    CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR();
    CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR(
        const VECTOR<T,D>& beta_dv_over_dx_dx,
        const unsigned char* p_value,
        BEGIN_TAG);
    CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR(
        const VECTOR<T,D>& beta_dv_over_dx_dx,
        const unsigned char* p_value,
        END_TAG);

    bool Valid() const;
    void Next();
    void Prev();

private:
    VECTOR<T,D> m_beta_dv_over_dx_dx;
    int m_index;
    const unsigned char* mp_value;

    friend class boost::iterator_core_access;

    reference dereference() const;
    bool equal(const CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR& other) const;
    void increment();
    void decrement();
    void advance(difference_type n);
    difference_type distance_to(const CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR& other) const;
};

//##################################################################### 
//##################################################################### 

template< class T, int D >
inline
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR()
{ }

template< class T, int D >
inline
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR(
    const VECTOR<T,D>& beta_dv_over_dx_dx,
    const unsigned char* p_value,
    BEGIN_TAG)
    : m_beta_dv_over_dx_dx(beta_dv_over_dx_dx),
      m_index(0),
      mp_value(p_value)
{ }

template< class T, int D >
inline
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR(
    const VECTOR<T,D>& beta_dv_over_dx_dx,
    const unsigned char* p_value,
    END_TAG)
    : m_beta_dv_over_dx_dx(beta_dv_over_dx_dx),
      m_index(2*D + 1),
      mp_value(p_value)
{ }

template< class T, int D >
inline bool
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::
Valid() const
{ return m_index <= 2*D; }

template< class T, int D >
inline void
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::
Next()
{ increment(); }

template< class T, int D >
inline void
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::
Prev()
{ decrement(); }

template< class T, int D >
inline typename CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::reference
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::
dereference() const
{ return (m_index == D ? m_beta_dv_over_dx_dx.Sum() : -m_beta_dv_over_dx_dx[D - std::abs(m_index - D) + 1]) * *mp_value; }

template< class T, int D >
inline bool
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::
equal(const CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR& other) const
{ return mp_value == other.mp_value; }

template< class T, int D >
inline void
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::
increment()
{
    assert(m_index <= 2*D);
    ++m_index;
    ++mp_value;
}

template< class T, int D >
inline void
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::
decrement()
{
    --m_index;
    --mp_value;
    assert(0 <= m_index);
}

template< class T, int D >
inline void
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::
advance(difference_type n)
{
    m_index += n;
    mp_value += n;
}

template< class T, int D >
inline typename CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::difference_type
CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR<T,D>::
distance_to(const CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR& other) const
{
    assert(other.m_index - m_index == other.mp_value - mp_value);
    return other.m_index - m_index;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_CROSS_CONSTBETA_STENCIL_VALUE_ITERATOR_HPP
