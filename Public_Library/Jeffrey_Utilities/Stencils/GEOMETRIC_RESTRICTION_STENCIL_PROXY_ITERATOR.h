//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR_HPP

#include <boost/iterator/iterator_traits.hpp>

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/ITERATOR_ADAPTOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE_ITERATOR.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int D >
class GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR;

template< class T, int D >
struct GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR_TRAITS
{
    typedef VECTOR<int,D> INDEX_TYPE;
    typedef T SCALAR_TYPE;
    typedef ITERATOR_ADAPTOR<
        GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>, // Derived
        MULTI_INDEX_CUBE_ITERATOR<D,-1,+1>,                // Base
        INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE >,            // Value
        boost::use_default,                                // CategoryOrTraversal
        INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE >,            // Reference
        boost::use_default                                 // Difference
    > ITERATOR_ADAPTOR_;
};

template< class T, int D >
class GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR
    : public GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR_TRAITS<T,D>::ITERATOR_ADAPTOR_
{
    typedef GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR_TRAITS<T,D> GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR_TRAITS_;
    typedef typename GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR_TRAITS_::ITERATOR_ADAPTOR_ ITERATOR_ADAPTOR_;
public:
    typedef typename ITERATOR_ADAPTOR_::value_type value_type;
    typedef typename ITERATOR_ADAPTOR_::reference reference;
    typedef typename ITERATOR_ADAPTOR_::difference_type difference_type;

    typedef typename ITERATOR_ADAPTOR_::base_type base_type;

    typedef typename GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR_TRAITS_::INDEX_TYPE INDEX_TYPE;
    typedef typename GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR_TRAITS_::SCALAR_TYPE SCALAR_TYPE;

    GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR();
    explicit GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& coarse_base_multi_index);
    GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& coarse_base_multi_index, BEGIN_TAG);
    GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& coarse_base_multi_index, END_TAG);

    INDEX_TYPE Fine_Base_Multi_Index() const;
    INDEX_TYPE Coarse_Base_Multi_Index() const;

    using ITERATOR_ADAPTOR_::base;

    INDEX_TYPE Index() const;
    SCALAR_TYPE Value() const;

    bool Valid() const;
    void Next();
    void Prev();

private:
    friend class boost::iterator_core_access;

    reference dereference() const;
};

//#####################################################################
//#####################################################################

template< class T, int D >
inline
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR()
{ }

template< class T, int D >
inline
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& coarse_base_multi_index)
    : ITERATOR_ADAPTOR_(base_type(2 * coarse_base_multi_index - 1))
{ }

template< class T, int D >
inline
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& coarse_base_multi_index, BEGIN_TAG)
    : ITERATOR_ADAPTOR_(base_type(2 * coarse_base_multi_index - 1, BEGIN_TAG()))
{ }

template< class T, int D >
inline
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR(const INDEX_TYPE& coarse_base_multi_index, END_TAG)
    : ITERATOR_ADAPTOR_(base_type(2 * coarse_base_multi_index - 1, END_TAG()))
{ }

template< class T, int D >
inline typename GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::INDEX_TYPE
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::
Fine_Base_Multi_Index() const
{ return base().Multi_Index_Cube().base_multi_index; }

template< class T, int D >
inline typename GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::INDEX_TYPE
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::
Coarse_Base_Multi_Index() const
{ return (Fine_Base_Multi_Index() + 1) / 2; }

template< class T, int D >
inline typename GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::INDEX_TYPE
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::
Index() const
{ return *base(); }

template< class T, int D >
inline typename GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::SCALAR_TYPE
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::
Value() const
{ return static_cast< SCALAR_TYPE >(1 << (D - base().Multi_Offset().L1_Norm())); }

template< class T, int D >
inline bool
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::
Valid() const
{ return base().Valid(); }

template< class T, int D >
inline void
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::
Next()
{ ITERATOR_ADAPTOR_::operator++(); }

template< class T, int D >
inline void
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::
Prev()
{ ITERATOR_ADAPTOR_::operator--(); }

template< class T, int D >
inline typename GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::reference
GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D>::
dereference() const
{ return reference(Index(), Value()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR_HPP
