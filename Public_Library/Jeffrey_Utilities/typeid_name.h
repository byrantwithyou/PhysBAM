//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_TYPEID_NAME_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_TYPEID_NAME_HPP

#include <typeinfo>

namespace PhysBAM
{

template< class T >
inline const char*
typeid_name()
{ return typeid(T).name(); }

template<>
inline const char*
typeid_name< float >()
{ return "float"; }

template<>
inline const char*
typeid_name< double >()
{ return "double"; }

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_TYPEID_NAME_HPP
