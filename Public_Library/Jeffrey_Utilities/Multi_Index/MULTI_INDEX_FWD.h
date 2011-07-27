//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_FWD_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_FWD_HPP

namespace PhysBAM
{

template< int D >
struct MULTI_INDEX_BOUND;

template< int D >
struct MULTI_INDEX_BOX;

template< int D, int MIN_OFFSET_, int MAX_OFFSET_ = -MIN_OFFSET_ >
struct MULTI_INDEX_CUBE;

template< int D, int MIN_INDEX_, int MAX_INDEX_ = -MIN_INDEX_ >
struct STATIC_MULTI_INDEX_CUBE;

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MULTI_INDEX_MULTI_INDEX_FWD_HPP
