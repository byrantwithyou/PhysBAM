//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_IDENTITY_STENCIL_PROXY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_IDENTITY_STENCIL_PROXY_HPP

#include <Jeffrey_Utilities/BEGIN_END_TAG.h>
#include <Jeffrey_Utilities/Stencils/IDENTITY_STENCIL_PROXY_ITERATOR.h>
#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>

namespace PhysBAM
{

template< class T_INDEX, class T >
struct IDENTITY_STENCIL_PROXY
{
    typedef T_INDEX INDEX_TYPE;
    typedef T SCALAR_TYPE;
    typedef INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE > INDEX_VALUE_TYPE;

    const INDEX_TYPE index;

    explicit IDENTITY_STENCIL_PROXY(const INDEX_TYPE& index_);

    static int N_Nonzero();

    typedef IDENTITY_STENCIL_PROXY_ITERATOR< T_INDEX, T > iterator;
    typedef iterator const_iterator;
    typedef typename iterator::reference reference;
    iterator begin() const;
    iterator end() const;
};

//#####################################################################
//#####################################################################

template< class T_INDEX, class T >
inline
IDENTITY_STENCIL_PROXY< T_INDEX, T >::
IDENTITY_STENCIL_PROXY(const INDEX_TYPE& index_)
    : index(index_)
{ }

template< class T_INDEX, class T >
inline int
IDENTITY_STENCIL_PROXY< T_INDEX, T >::
N_Nonzero()
{ return 1; }

template< class T_INDEX, class T >
inline typename IDENTITY_STENCIL_PROXY< T_INDEX, T >::iterator
IDENTITY_STENCIL_PROXY< T_INDEX, T >::
begin() const
{ return iterator(index, BEGIN_TAG()); }

template< class T_INDEX, class T >
inline typename IDENTITY_STENCIL_PROXY< T_INDEX, T >::iterator
IDENTITY_STENCIL_PROXY< T_INDEX, T >::
end() const
{ return iterator(index, END_TAG()); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_IDENTITY_STENCIL_PROXY_HPP
