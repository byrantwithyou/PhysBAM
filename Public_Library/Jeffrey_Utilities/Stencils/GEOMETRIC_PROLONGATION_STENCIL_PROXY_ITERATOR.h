//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR_HPP

#include <cassert>

#include <boost/iterator/iterator_traits.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/ITERATOR_ADAPTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOX.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOX_ITERATOR.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>

namespace PhysBAM
{

template< class T, int D >
class GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR;

template< class T, int D >
struct GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR_TRAITS
{
    typedef VECTOR<int,D> INDEX_TYPE;
    typedef T SCALAR_TYPE;
    typedef ITERATOR_ADAPTOR<
        GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>, // Derived
        MULTI_INDEX_BOX_ITERATOR<D>,                        // Base
        INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE >,             // Value
        boost::use_default,                                 // CategoryOrTraversal
        INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE >,             // Reference
        boost::use_default                                  // Difference
    > ITERATOR_ADAPTOR_;
};

template< class T, int D >
class GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR
    : public GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR_TRAITS<T,D>::ITERATOR_ADAPTOR_
{
    typedef GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR_TRAITS<T,D> GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR_TRAITS_;
    typedef typename GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR_TRAITS_::ITERATOR_ADAPTOR_ ITERATOR_ADAPTOR_;
public:
    typedef typename ITERATOR_ADAPTOR_::value_type value_type;
    typedef typename ITERATOR_ADAPTOR_::reference reference;
    typedef typename ITERATOR_ADAPTOR_::difference_type difference_type;

    typedef typename ITERATOR_ADAPTOR_::base_type base_type;

    typedef typename GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR_TRAITS_::INDEX_TYPE INDEX_TYPE;
    typedef typename GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR_TRAITS_::SCALAR_TYPE SCALAR_TYPE;

    typedef MULTI_INDEX_BOX<D> MULTI_INDEX_BOX_TYPE;

    GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR();
    explicit GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& fine_base_multi_index);
    GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& fine_base_multi_index, BEGIN_TAG);
    GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& fine_base_multi_index, END_TAG);

    using ITERATOR_ADAPTOR_::base;

    INDEX_TYPE Index() const;
    SCALAR_TYPE Value() const;

    bool Valid() const;
    void Next();
    void Prev();

private:
    const T m_value;

    static MULTI_INDEX_BOX_TYPE Coarse_Multi_Index_Box(const INDEX_TYPE fine_base_multi_index);
    static T Value(const INDEX_TYPE fine_base_multi_index);

    friend class boost::iterator_core_access;

    reference dereference() const;
};

//#####################################################################
//#####################################################################

template< class T, int D >
inline
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR()
{ }

template< class T, int D >
inline
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& fine_base_multi_index)
    : ITERATOR_ADAPTOR_(base_type(Coarse_Multi_Index_Box(fine_base_multi_index))),
      m_value(Value(fine_base_multi_index))
{ }

template< class T, int D >
inline
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& fine_base_multi_index, BEGIN_TAG)
    : ITERATOR_ADAPTOR_(base_type(Coarse_Multi_Index_Box(fine_base_multi_index), BEGIN_TAG())),
      m_value(Value(fine_base_multi_index))
{ }

template< class T, int D >
inline
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& fine_base_multi_index, END_TAG)
    : ITERATOR_ADAPTOR_(base_type(Coarse_Multi_Index_Box(fine_base_multi_index), END_TAG())),
      m_value(Value(fine_base_multi_index))
{ }

template< class T, int D >
inline typename GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::INDEX_TYPE
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::
Index() const
{ return *base(); }

template< class T, int D >
inline typename GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::SCALAR_TYPE
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::
Value() const
{ return m_value; }

template< class T, int D >
inline bool
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::
Valid() const
{ return base().Valid(); }

template< class T, int D >
inline void
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::
Next()
{ ITERATOR_ADAPTOR_::operator++(); }

template< class T, int D >
inline void
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::
Prev()
{ ITERATOR_ADAPTOR_::operator--(); }

template< class T, int D >
inline typename GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::MULTI_INDEX_BOX_TYPE
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::
Coarse_Multi_Index_Box(const INDEX_TYPE fine_base_multi_index)
{
    MULTI_INDEX_BOX_TYPE coarse_multi_index_box((fine_base_multi_index + 1) / 2);
    for(int d = 1; d <= D; ++d)
        if(fine_base_multi_index[d] & 1 == 0)
            ++coarse_multi_index_box.max_multi_index[d];
    return coarse_multi_index_box;
}

template< class T, int D >
inline T
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::
Value(const INDEX_TYPE fine_base_multi_index)
{
    int value_inv = 1;
    for(int d = 1; d <= D; ++d)
        if(fine_base_multi_index[d] & 1 == 0)
            value_inv *= 2;
    return static_cast<T>(1) / value_inv;
}

template< class T, int D >
inline typename GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::reference
GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D>::
dereference() const
{ return reference(Index(), Value()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR_HPP
