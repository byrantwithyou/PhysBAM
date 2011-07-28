//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_RESTRICTION_STENCIL_PROXY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_RESTRICTION_STENCIL_PROXY_HPP

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Stencils/GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int D >
struct GEOMETRIC_RESTRICTION_STENCIL_PROXY
{
    typedef VECTOR<int,D> INDEX_TYPE;
    typedef T SCALAR_TYPE;
    typedef void STENCIL_TYPE;

    static const int DIMENSION = D;
    static const int WIDTH = 3;

    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS(
        GEOMETRIC_RESTRICTION_STENCIL_PROXY,
        (( typename INDEX_TYPE const, coarse_base_multi_index ))
    )

    static int N_Nonzero();

    typedef GEOMETRIC_RESTRICTION_STENCIL_PROXY_ITERATOR<T,D> iterator;
    typedef iterator const_iterator;
    typedef typename iterator::reference reference;
    iterator begin() const;
    iterator end() const;
};

//#####################################################################
//#####################################################################

template< class T, int D >
inline int
GEOMETRIC_RESTRICTION_STENCIL_PROXY<T,D>::
N_Nonzero()
{ return STATIC_POW_C< 3, DIMENSION >::value; }

template< class T, int D >
inline typename GEOMETRIC_RESTRICTION_STENCIL_PROXY<T,D>::iterator
GEOMETRIC_RESTRICTION_STENCIL_PROXY<T,D>::
begin() const
{ return iterator(coarse_base_multi_index, BEGIN_TAG()); }

template< class T, int D >
inline typename GEOMETRIC_RESTRICTION_STENCIL_PROXY<T,D>::iterator
GEOMETRIC_RESTRICTION_STENCIL_PROXY<T,D>::
end() const
{ return iterator(coarse_base_multi_index, END_TAG()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_RESTRICTION_STENCIL_PROXY_HPP
