//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_PROLONGATION_STENCIL_PROXY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_PROLONGATION_STENCIL_PROXY_HPP

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Stencils/GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int D >
struct GEOMETRIC_PROLONGATION_STENCIL_PROXY
{
    typedef VECTOR<int,D> INDEX_TYPE;
    typedef T SCALAR_TYPE;
    typedef void STENCIL_TYPE;

    static const int DIMENSION = D;

    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_MEMBERS(
        GEOMETRIC_PROLONGATION_STENCIL_PROXY,
        (( typename INDEX_TYPE const, fine_base_multi_index ))
    )

    int N_Nonzero() const;

    typedef GEOMETRIC_PROLONGATION_STENCIL_PROXY_ITERATOR<T,D> iterator;
    typedef iterator const_iterator;
    typedef typename iterator::reference reference;
    iterator begin() const;
    iterator end() const;
};

//#####################################################################
//#####################################################################

template< class T, int D >
inline int
GEOMETRIC_PROLONGATION_STENCIL_PROXY<T,D>::
N_Nonzero() const
{
    int n_nonzero = 1;
    for(int d = 1; d <= D; ++d)
        if(fine_base_multi_index[d] & 1 == 0)
            n_nonzero *= 2;
    return n_nonzero;
}

template< class T, int D >
inline typename GEOMETRIC_PROLONGATION_STENCIL_PROXY<T,D>::iterator
GEOMETRIC_PROLONGATION_STENCIL_PROXY<T,D>::
begin() const
{ return iterator(fine_base_multi_index, BEGIN_TAG()); }

template< class T, int D >
inline typename GEOMETRIC_PROLONGATION_STENCIL_PROXY<T,D>::iterator
GEOMETRIC_PROLONGATION_STENCIL_PROXY<T,D>::
end() const
{ return iterator(fine_base_multi_index, END_TAG()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_GEOMETRIC_PROLONGATION_STENCIL_PROXY_HPP
