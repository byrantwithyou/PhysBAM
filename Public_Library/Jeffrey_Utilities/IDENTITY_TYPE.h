//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_IDENTITY_TYPE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_IDENTITY_TYPE_HPP

#define PHYSBAM_IDENTITY_TYPE( ParenthesizedT ) \
    PhysBAM::Detail_IDENTITY_TYPE::IDENTITY_TYPE_IMPL< void ParenthesizedT >::type

namespace PhysBAM
{

namespace Detail_IDENTITY_TYPE
{

template< class >
struct IDENTITY_TYPE_IMPL;
template< class T >
struct IDENTITY_TYPE_IMPL< void ( T ) >
{ typedef T type; };

} // namespace Detail_IDENTITY_TYPE

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_IDENTITY_TYPE_HPP
