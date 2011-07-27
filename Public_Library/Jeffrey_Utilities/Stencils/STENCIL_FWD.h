//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_STENCIL_FWD_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_STENCIL_FWD_HPP

namespace PhysBAM
{

template< int D >
struct CROSS_CONSTBETA_STENCIL;

template< class T, int D >
struct CROSS_STENCIL;

template< class T, int D, int MIN_OFFSET_, int MAX_OFFSET_ = -MIN_OFFSET_ >
struct CUBE_STENCIL;

template< class T_INDEX, class T, int MAX_N_NONZERO = -1 >
struct UNSTRUCTURED_STENCIL;

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_STENCIL_FWD_HPP
