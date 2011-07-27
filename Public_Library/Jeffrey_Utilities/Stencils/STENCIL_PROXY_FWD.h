//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_STENCIL_PROXY_FWD_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_STENCIL_PROXY_FWD_HPP

namespace PhysBAM
{

template< class T_STENCIL, class T >
struct CROSS_CONSTBETA_STENCIL_PROXY;

template< class T_STENCIL >
struct CROSS_STENCIL_PROXY;

template< class T_STENCIL >
struct CUBE_STENCIL_PROXY;

template< class T, int D >
struct GEOMETRIC_PROLONGATION_STENCIL_PROXY;

template< class T, int D >
struct GEOMETRIC_RESTRICTION_STENCIL_PROXY;

template<
    class T_BASE,
    class T_INDEX_TRANSFORM,
    class T_INVERSE_INDEX_TRANSFORM = T_INDEX_TRANSFORM
>
struct INDEX_TRANSFORM_STENCIL_PROXY;

template< class T_BASE >
struct SCALED_STENCIL_PROXY;

// TODO: REMOVE
template< class T_BASE, class T_INDEX, class T_INDEX_CONVERTER >
struct STENCIL_INDEX_ADAPTOR_PROXY;

template< class T_STENCIL >
struct UNSTRUCTURED_STENCIL_PROXY;

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_STENCILS_STENCIL_PROXY_FWD_HPP
