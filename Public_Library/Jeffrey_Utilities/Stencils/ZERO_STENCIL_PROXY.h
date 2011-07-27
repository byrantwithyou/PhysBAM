//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_ZERO_STENCIL_PROXY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_ZERO_STENCIL_PROXY_HPP

#include <Jeffrey_Utilities/Stencils/INDEX_VALUE.h>

namespace PhysBAM
{

template< class T_INDEX, class T >
struct ZERO_STENCIL_PROXY
{
    typedef T_INDEX INDEX_TYPE;
    typedef T SCALAR_TYPE;
    typedef INDEX_VALUE< INDEX_TYPE, SCALAR_TYPE > INDEX_VALUE_TYPE;

    static int N_Nonzero();

    typedef const INDEX_VALUE_TYPE* iterator;
    typedef iterator const_iterator;
    typedef const INDEX_VALUE_TYPE& reference;
    iterator begin() const;
    iterator end() const;
};

//#####################################################################
//#####################################################################

template< class T_INDEX, class T >
inline int
ZERO_STENCIL_PROXY< T_INDEX, T >::
N_Nonzero()
{ return 0; }

template< class T_INDEX, class T >
inline typename ZERO_STENCIL_PROXY< T_INDEX, T >::iterator
ZERO_STENCIL_PROXY< T_INDEX, T >::
begin() const
{ return static_cast< iterator >(0); }

template< class T_INDEX, class T >
inline typename ZERO_STENCIL_PROXY< T_INDEX, T >::iterator
ZERO_STENCIL_PROXY< T_INDEX, T >::
end() const
{ return static_cast< iterator >(0); }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_ZERO_STENCIL_PROXY_HPP
